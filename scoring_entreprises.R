# =============================================================================
#  scoring_entreprises.R — modèle de scoring du risque de défaillance des entreprises
# =============================================================================
#  Cible     : défaillance (procédure collective, défaut bancaire > 90 jours…)
#              dans les 12 mois suivant l'arrêté des comptes.
#  Entrées   : liasse fiscale simplifiée (bilan + compte de résultat), secteur,
#              âge, effectif, incidents de paiement.
#  Sorties   : PD à 1 an, score en points, classe de risque 1 (meilleure) à 8,
#              trois motifs explicatifs par entreprise.
#
#  Lancer    : Rscript scoring_entreprises.R
#  Vos données : renseigner FICHIER_DONNEES (CSV ; séparateur « ; », décimale « , »)
#              avec les colonnes décrites dans `simuler_entreprises()`.
#              Sinon un jeu réaliste de 15 000 bilans (2019-2023) est simulé.
# =============================================================================

FICHIER_DONNEES <- NULL          # ex. "mes_bilans.csv"
DOSSIER_SORTIE  <- "sorties/entreprises"
ANNEE_HORS_PERIODE <- 2023       # validation out-of-time : dernière année
TAUX_DEFAUT_LONG_TERME <- NA     # ex. 0.045 : recale les PD sur une moyenne de cycle (NA = pas de recalage)
set.seed(2026)

dossier_script <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
racine <- dossier_script()
source(file.path(racine, "fonctions_scoring.R"))
sortie <- file.path(racine, DOSSIER_SORTIE)
dir.create(sortie, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 1. Données
# -----------------------------------------------------------------------------

#' Jeu simulé : bilans et comptes de résultat cohérents, avec valeurs manquantes,
#' valeurs extrêmes, effet sectoriel, choc 2020 (aides publiques) et hausse des
#' défaillances en 2023 (dérive à détecter en hors période).
simuler_entreprises <- function(n = 15000) {
  secteurs <- c("Industrie", "Commerce", "BTP", "Services", "Transport", "Hôtellerie-restauration")
  d <- data.frame(
    id      = sprintf("E%05d", seq_len(n)),
    annee   = sample(2019:2023, n, replace = TRUE),
    secteur = sample(secteurs, n, replace = TRUE, prob = c(.18, .24, .16, .25, .07, .10)),
    age     = round(rexp(n, 1 / 14) + runif(n, 0, 2), 1),
    effectif = pmax(1, round(exp(rnorm(n, 2.3, 1.2))))
  )
  z <- rnorm(n)                                   # santé « latente » de l'entreprise
  rot <- c(Industrie = .9, Commerce = .55, BTP = .8, Services = .7, Transport = 1, `Hôtellerie-restauration` = 1.2)[d$secteur]
  pmarge <- c(Industrie = .08, Commerce = .04, BTP = .06, Services = .10, Transport = .07, `Hôtellerie-restauration` = .09)[d$secteur]
  d$chiffre_affaires <- round(d$effectif * 120 * exp(rnorm(n, 0, .5)))                    # k€
  d$ebe            <- round(d$chiffre_affaires * (pmarge + .035 * z + rnorm(n, 0, .05)))
  d$total_actif    <- round(d$chiffre_affaires * rot * exp(rnorm(n, 0, .3)))
  d$fonds_propres  <- round(d$total_actif * (.30 + .10 * z + rnorm(n, 0, .13)))
  d$dettes_financieres <- round(d$total_actif * pmax(0, .32 - .07 * z + rnorm(n, 0, .14)))
  d$tresorerie     <- round(d$total_actif * pmax(0, .08 + .035 * z + rnorm(n, 0, .06)))
  d$resultat_net   <- round(d$ebe - .05 * d$total_actif - .035 * d$dettes_financieres)
  d$actif_circulant <- round(d$total_actif * runif(n, .35, .7))
  d$passif_circulant <- round(d$actif_circulant / exp(.25 + .22 * z + rnorm(n, 0, .25)))
  dso <- c(Industrie = 60, Commerce = 25, BTP = 75, Services = 55, Transport = 50, `Hôtellerie-restauration` = 8)[d$secteur]
  d$creances_clients <- round(d$chiffre_affaires * pmax(0, dso - 9 * z + rnorm(n, 0, 12)) / 365)
  d$ca_precedent   <- round(d$chiffre_affaires / exp(.03 + .06 * z - .12 * (d$annee == 2020) + rnorm(n, 0, .10)))
  d$incidents_paiement <- rpois(n, exp(-2.2 - .9 * z))
  # valeurs manquantes : comptes partiels, entreprises jeunes sans exercice précédent
  d$ca_precedent[d$age < 1 | runif(n) < .03] <- NA
  d$creances_clients[runif(n) < .06] <- NA
  d$passif_circulant[runif(n) < .03] <- NA
  # valeurs aberrantes (erreurs de saisie)
  i <- sample.int(n, 30); d$chiffre_affaires[i] <- d$chiffre_affaires[i] * 1000

  # processus de défaut (inconnu du modélisateur)
  eta <- -3.75 - 1.05 * z + .45 * pmin(d$incidents_paiement, 4) + .7 * (d$fonds_propres < 0) +
    .5 * (d$age < 3) - .15 * log(d$effectif) +
    c(Industrie = 0, Commerce = .1, BTP = .35, Services = -.1, Transport = .15, `Hôtellerie-restauration` = .45)[d$secteur] +
    c(`2019` = 0, `2020` = -.45, `2021` = -.25, `2022` = .1, `2023` = .35)[as.character(d$annee)]
  d$defaut <- rbinom(n, 1, plogis(eta))
  d
}

#' Ratios financiers à partir des agrégats comptables. Les divisions par zéro
#' et les dénominateurs négatifs sont traités explicitement.
calculer_ratios <- function(d) {
  div <- function(a, b) ifelse(is.finite(a / b) & b != 0, a / b, NA)
  data.frame(
    id = d$id, annee = d$annee,
    defaut = if (is.null(d$defaut)) NA else d$defaut,   # inconnu pour une nouvelle entreprise
    secteur = d$secteur,
    age = d$age,
    log_effectif = log(d$effectif),
    autonomie_financiere = div(d$fonds_propres, d$total_actif),
    # capacité de remboursement (années d'EBE) : EBE <= 0 -> plafond, le plus risqué
    capacite_remboursement = ifelse(d$ebe > 0, pmin(div(d$dettes_financieres - d$tresorerie, d$ebe), 30), 30),
    marge_ebe = div(d$ebe, d$chiffre_affaires),
    rentabilite_actif = div(d$resultat_net, d$total_actif),
    liquidite_generale = div(d$actif_circulant, d$passif_circulant),
    tresorerie_actif = div(d$tresorerie, d$total_actif),
    delai_clients_jours = div(d$creances_clients, d$chiffre_affaires) * 365,
    croissance_ca = div(d$chiffre_affaires, d$ca_precedent) - 1,
    rotation_actif = div(d$chiffre_affaires, d$total_actif),
    incidents_paiement = d$incidents_paiement
  )
}

brut <- if (is.null(FICHIER_DONNEES)) simuler_entreprises() else read.csv2(FICHIER_DONNEES, stringsAsFactors = FALSE)
donnees <- calculer_ratios(brut)
candidates <- setdiff(names(donnees), c("id", "annee", "defaut"))

message(sprintf("%d bilans, taux de défaut %.2f %%", nrow(donnees), 100 * mean(donnees$defaut)))
print(round(100 * tapply(donnees$defaut, donnees$annee, mean), 2))


# -----------------------------------------------------------------------------
# 2. Échantillons : apprentissage / test (même période) / hors période
# -----------------------------------------------------------------------------
periode <- donnees[donnees$annee <  ANNEE_HORS_PERIODE, ]
hors_periode <- donnees[donnees$annee >= ANNEE_HORS_PERIODE, ]
i_app <- decoupe_stratifiee(periode$defaut, 0.7)
app <- periode[i_app, ]; test <- periode[-i_app, ]
message(sprintf("Apprentissage %d | test %d | hors période %d", nrow(app), nrow(test), nrow(hors_periode)))


# -----------------------------------------------------------------------------
# 3. Modèle
# -----------------------------------------------------------------------------
parametres <- list(part_min = 0.05, nb_max = 7, iv_min = 0.02, cor_max = 0.6, vif_max = 5, p_max = 0.05,
                   score_ref = 600, cote_ref = 50, pdo = 20)
message("Ajustement de la grille :")
modele <- ajuster_score(app, "defaut", candidates, parametres)
print(modele)

# Le taux de défaut de l'apprentissage (2019-2022) intègre les années d'aides
# publiques : il sous-estime le risque d'une année normale. Pour une PD
# « à travers le cycle », on recale la constante sur un taux de long terme.
if (!is.na(TAUX_DEFAUT_LONG_TERME)) modele <- recalibrer(modele, TAUX_DEFAUT_LONG_TERME, app)

# échelle de notation à taux de défaut monotone
score_app <- predire_score(modele, app, 1)$score
modele$classes <- construire_classes(score_app, app$defaut, nb_classes = 8)
if (!is.na(TAUX_DEFAUT_LONG_TERME))
  modele$classes$pd <- tapply(predire_score(modele, app, 1)$pd, attribuer_classe(modele$classes, score_app), mean)[as.character(modele$classes$classe)]
print(modele$classes, row.names = FALSE)


# -----------------------------------------------------------------------------
# 4. Validation
# -----------------------------------------------------------------------------
pred <- lapply(list(Apprentissage = app, Test = test, `Hors période` = hors_periode),
               function(d) cbind(predire_score(modele, d), y = d$defaut))

perf <- t(sapply(pred, function(p) performances(p$pd, p$y)))
message("\nPerformances :"); print(round(perf, 4))

boot <- t(sapply(pred[-1], function(p) gini_bootstrap(p$pd, p$y, B = 500)))
message("Gini avec IC bootstrap à 95 % :"); print(round(100 * boot, 1))

bt <- backtest_classes(modele$classes, pred$`Hors période`$classe, pred$`Hors période`$y)
message("Backtesting par classe (hors période) :"); print(bt, digits = 3, row.names = FALSE)

alerte_calibration(perf)          # le cas échéant : renseigner TAUX_DEFAUT_LONG_TERME

stab <- stabilite(modele, app, hors_periode)
message("Stabilité apprentissage -> hors période (PSI) :"); print(stab, digits = 3, row.names = FALSE)

message("Validation croisée répétée (3 x 5 plis, toute la chaîne réapprise)…")
vc <- validation_croisee(periode, "defaut", candidates, parametres, k = 5, repetitions = 3)
message(sprintf("  Gini = %.1f %% ± %.1f", 100 * vc$gini_moyen, 100 * vc$gini_ecart_type))
print(vc$frequence_selection, digits = 2, row.names = FALSE)

message("Modèles challengers :")
rf <- challenger_foret(app, list(test = test, hp = hors_periode), "defaut", candidates)
en <- challenger_glmnet(modele, app, list(test = test, hp = hors_periode))
challengers <- data.frame(
  modele = c("Grille de score (logistique WoE)", "Elastic net (WoE)", "Forêt aléatoire (brut)"),
  gini_test = c(gini(pred$Test$pd, test$defaut),
                if (!is.null(en)) gini(en$pd$test, test$defaut) else NA,
                if (!is.null(rf)) gini(rf$pd$test, test$defaut) else NA),
  gini_hors_periode = c(gini(pred$`Hors période`$pd, hors_periode$defaut),
                        if (!is.null(en)) gini(en$pd$hp, hors_periode$defaut) else NA,
                        if (!is.null(rf)) gini(rf$pd$hp, hors_periode$defaut) else NA))
print(challengers, digits = 3, row.names = FALSE)


# -----------------------------------------------------------------------------
# 5. Exports
# -----------------------------------------------------------------------------
write.csv2(modele$grille, file.path(sortie, "grille_score.csv"), row.names = FALSE)
write.csv2(do.call(rbind, lapply(modele$bins, `[[`, "table")), file.path(sortie, "binning_toutes_variables.csv"), row.names = FALSE)
write.csv2(modele$classes, file.path(sortie, "classes_de_risque.csv"), row.names = FALSE)
write.csv2(data.frame(echantillon = rownames(perf), perf), file.path(sortie, "performances.csv"), row.names = FALSE)
write.csv2(bt, file.path(sortie, "backtesting_hors_periode.csv"), row.names = FALSE)
write.csv2(stab, file.path(sortie, "stabilite_psi.csv"), row.names = FALSE)
write.csv2(vc$frequence_selection, file.path(sortie, "validation_croisee_selection.csv"), row.names = FALSE)
write.csv2(challengers, file.path(sortie, "challengers.csv"), row.names = FALSE)
write.csv2(cbind(hors_periode[c("id", "annee", "secteur")], pred$`Hors période`[setdiff(names(pred$`Hors période`), "y")]),
           file.path(sortie, "scores_hors_periode.csv"), row.names = FALSE)
graphiques(file.path(sortie, "graphiques.pdf"), modele,
           lapply(pred, function(p) list(pd = p$pd, y = p$y, score = p$score)), vc)
saveRDS(modele, file.path(sortie, "modele_entreprises.rds"))
message("\nRésultats écrits dans ", sortie)

# Utilisation en production :
#   source("fonctions_scoring.R")
#   modele <- readRDS("sorties/entreprises/modele_entreprises.rds")
#   predire_score(modele, calculer_ratios(nouveaux_bilans))
