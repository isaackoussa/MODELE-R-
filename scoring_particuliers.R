# =============================================================================
#  scoring_particuliers.R — score d'octroi de crédit à la consommation (particuliers)
# =============================================================================
#  Cible     : défaut à 12 mois (3 échéances impayées ou plus, déchéance du terme…).
#  Entrées   : situation professionnelle, revenus et charges, logement, relation
#              bancaire, comportement de paiement, caractéristiques du crédit.
#  Sorties   : PD, score en points, classe de risque, décision (seuil optimisé
#              sur la rentabilité), trois motifs de refus par demande.
#
#  Cadre     : le scoring de crédit des personnes physiques est un système « à haut
#              risque » (AI Act, annexe III) et une décision automatisée encadrée
#              (RGPD art. 22) : le modèle doit être explicable, ne pas utiliser de
#              critère discriminatoire et faire l'objet d'un audit d'équité. D'où :
#                - attributs protégés EXCLUS du modèle, conservés uniquement pour l'audit ;
#                - motifs de refus individuels ;
#                - audit d'équité (ratio d'impact, calibration et AUC par groupe) ;
#                - recherche de variables proxy des attributs protégés.
#
#  Lancer    : Rscript scoring_particuliers.R
#  Vos données : renseigner FICHIER_DONNEES (CSV ; « ; » et décimale « , ») avec
#              les colonnes de `simuler_particuliers()`. Sinon 20 000 dossiers simulés.
# =============================================================================

FICHIER_DONNEES <- NULL
DOSSIER_SORTIE  <- "sorties/particuliers"
# Critères de discrimination prohibés (C. pénal art. 225-1) : jamais dans le modèle.
# L'âge est un critère protégé lui aussi, admis en crédit s'il est objectivement
# justifié : ajoutez "age" ici pour l'exclure.
ATTRIBUTS_PROTEGES <- c("sexe", "nationalite", "situation_familiale")
DEBUT_HORS_PERIODE <- as.Date("2024-07-01")
MARGE_NETTE <- 0.06              # marge sur la durée de vie d'un bon dossier (en % du montant)
LGD         <- 0.55              # perte en cas de défaut (en % du montant)
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

#' Jeu simulé réaliste. Le défaut ne dépend PAS directement du sexe ni de la
#' nationalité, mais ceux-ci sont corrélés au revenu et au contrat : l'audit
#' doit donc rester vigilant aux effets indirects.
simuler_particuliers <- function(n = 20000) {
  d <- data.frame(id = sprintf("P%06d", seq_len(n)),
                  date_demande = as.Date("2022-01-01") + sort(sample.int(1095, n, replace = TRUE)))
  d$sexe <- sample(c("F", "H"), n, replace = TRUE)
  d$nationalite <- sample(c("Française", "UE", "Hors UE"), n, replace = TRUE, prob = c(.86, .07, .07))
  d$age <- round(pmin(85, 18 + rgamma(n, 3.2, 1 / 9)))
  d$situation_familiale <- ifelse(d$age < 26, sample(c("Célibataire", "Couple"), n, TRUE, c(.75, .25)),
                                  sample(c("Célibataire", "Couple", "Divorcé(e)/séparé(e)", "Veuf(ve)"), n, TRUE, c(.3, .52, .14, .04)))
  d$type_contrat <- ifelse(d$age >= 63 & runif(n) < .85, "Retraité",
                     sample(c("CDI", "Fonctionnaire", "CDD/intérim", "Indépendant", "Sans emploi"), n, TRUE,
                            c(.56, .14, .14, .10, .06)))
  z <- rnorm(n)                                          # solidité financière latente
  base_rev <- c(CDI = 2300, Fonctionnaire = 2350, `CDD/intérim` = 1650, Indépendant = 2500,
                `Sans emploi` = 1050, Retraité = 1800)[d$type_contrat]
  d$revenu_mensuel <- round(base_rev * exp(.35 * z + rnorm(n, 0, .3)) * ifelse(d$sexe == "F", .9, 1))
  d$anciennete_emploi_mois <- ifelse(d$type_contrat %in% c("Sans emploi", "Retraité"), NA,
                                     round(pmin(rexp(n, 1 / 70), 12 * (d$age - 17))))
  d$statut_logement <- ifelse(d$age < 25 & runif(n) < .4, "Hébergé",
                         sample(c("Propriétaire", "Accédant", "Locataire"), n, TRUE,
                                prob = c(.28, .30, .42)))
  d$loyer <- ifelse(d$statut_logement == "Locataire", round(d$revenu_mensuel * runif(n, .2, .4)),
             ifelse(d$statut_logement == "Accédant", 0, 0))
  d$charges_credits <- round(d$revenu_mensuel * pmax(0, rnorm(n, .15 - .04 * z, .10)) +
                             ifelse(d$statut_logement == "Accédant", d$revenu_mensuel * .25, 0))
  d$nb_credits_en_cours <- rpois(n, exp(.1 - .25 * z))
  d$montant_demande <- round(pmax(500, exp(rnorm(n, 9, .7))), -2)
  d$duree_mois <- sample(c(12, 24, 36, 48, 60, 72, 84), n, TRUE)
  mensualite <- d$montant_demande * (.07 / 12) / (1 - (1 + .07 / 12)^-d$duree_mois)
  d$anciennete_banque_mois <- round(pmin(rexp(n, 1 / 90), 12 * (d$age - 17)))
  d$jours_decouvert_6m <- pmin(180, rpois(n, exp(1.2 - .9 * z)))
  d$incidents_12m <- rpois(n, exp(-2.3 - .8 * z))
  d$utilisation_revolving <- ifelse(runif(n) < .35, pmin(1.2, plogis(-.5 - .9 * z + rnorm(n, 0, .8))), NA)
  d$epargne <- round(exp(rnorm(n, 7.8 + .9 * z, 1.4)) * (runif(n) > .15))

  # dérive 2024 : inflation -> charges plus lourdes
  infl <- d$date_demande >= as.Date("2024-01-01")
  d$charges_credits[infl] <- round(d$charges_credits[infl] * 1.08)

  # variables dérivées (calculées aussi en production)
  d$taux_endettement <- round((d$charges_credits + mensualite) / d$revenu_mensuel, 4)
  d$reste_a_vivre <- round(d$revenu_mensuel - d$charges_credits - mensualite - d$loyer)
  d$ratio_montant_revenu <- round(d$montant_demande / (12 * d$revenu_mensuel), 4)

  # processus de défaut (inconnu du modélisateur) — sans effet direct du sexe ni de la nationalité
  eta <- -3.3 - .75 * z + 1.6 * pmin(d$taux_endettement - .33, .4) * (d$taux_endettement > .33) +
    .35 * pmin(d$incidents_12m, 3) + .012 * d$jours_decouvert_6m - .002 * pmin(d$anciennete_banque_mois, 240) +
    c(CDI = 0, Fonctionnaire = -.4, `CDD/intérim` = .45, Indépendant = .35, `Sans emploi` = .9, Retraité = -.1)[d$type_contrat] +
    c(Propriétaire = -.35, Accédant = -.1, Locataire = .15, Hébergé = .2)[d$statut_logement] +
    .5 * (d$age < 25) + .2 * (d$duree_mois >= 72) + .35 * infl
  d$defaut <- rbinom(n, 1, plogis(eta))
  d$charges_credits <- d$loyer <- NULL
  d
}

donnees <- if (is.null(FICHIER_DONNEES)) simuler_particuliers() else {
  x <- read.csv2(FICHIER_DONNEES, stringsAsFactors = FALSE); x$date_demande <- as.Date(x$date_demande); x
}
candidates <- setdiff(names(donnees), c("id", "date_demande", "defaut", ATTRIBUTS_PROTEGES))
message(sprintf("%d dossiers, taux de défaut %.2f %%", nrow(donnees), 100 * mean(donnees$defaut)))
message("Variables candidates : ", paste(candidates, collapse = ", "))
message("Exclues (attributs protégés, audit seulement) : ", paste(ATTRIBUTS_PROTEGES, collapse = ", "))


# -----------------------------------------------------------------------------
# 2. Échantillons
# -----------------------------------------------------------------------------
periode <- donnees[donnees$date_demande < DEBUT_HORS_PERIODE, ]
hors_periode <- donnees[donnees$date_demande >= DEBUT_HORS_PERIODE, ]
i_app <- decoupe_stratifiee(periode$defaut, 0.7)
app <- periode[i_app, ]; test <- periode[-i_app, ]
message(sprintf("Apprentissage %d | test %d | hors période %d", nrow(app), nrow(test), nrow(hors_periode)))


# -----------------------------------------------------------------------------
# 3. Modèle
# -----------------------------------------------------------------------------
parametres <- list(part_min = 0.05, nb_max = 7, iv_min = 0.02, cor_max = 0.6, vif_max = 5, p_max = 0.05,
                   score_ref = 600, cote_ref = 30, pdo = 20)
message("Ajustement de la grille :")
modele <- ajuster_score(app, "defaut", candidates, parametres)
print(modele)
modele$classes <- construire_classes(predire_score(modele, app, 1)$score, app$defaut, nb_classes = 10)
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
alerte_calibration(perf)          # le cas échéant : modele <- recalibrer(modele, taux_cible, app)
stab <- stabilite(modele, app, hors_periode)
message("Stabilité (PSI) :"); print(stab, digits = 3, row.names = FALSE)

message("Validation croisée répétée (3 x 5 plis)…")
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
# 5. Politique d'octroi : seuil qui maximise la rentabilité espérée sur le test
# -----------------------------------------------------------------------------
strategie <- table_strategie(pred$Test$score, pred$Test$pd, test$montant_demande,
                             marge = MARGE_NETTE, lgd = LGD, y = test$defaut)
seuil <- strategie$seuil[which.max(strategie$rentabilite_esperee)]
message(sprintf("\nSeuil d'acceptation retenu : %d points", seuil))
print(strategie[abs(strategie$seuil - seuil) <= 15, ], digits = 3, row.names = FALSE)
modele$seuil_acceptation <- seuil


# -----------------------------------------------------------------------------
# 6. Audit d'équité et proxys (sur test + hors période)
# -----------------------------------------------------------------------------
audit_donnees <- rbind(test, hors_periode)
audit_pred <- rbind(pred$Test, pred$`Hors période`)
audit_donnees$tranche_age <- cut(audit_donnees$age, c(17, 25, 35, 50, 65, Inf),
                                 labels = c("18-25", "26-35", "36-50", "51-65", "66+"))
accepte <- audit_pred$score >= seuil
equite <- do.call(rbind, lapply(c(ATTRIBUTS_PROTEGES, "tranche_age"), function(a)
  cbind(attribut = a, audit_equite(audit_donnees[[a]], accepte, audit_pred$y, audit_pred$pd))))
message("\nAudit d'équité :"); print(equite, digits = 3, row.names = FALSE)
proxys <- do.call(rbind, lapply(ATTRIBUTS_PROTEGES, function(a) detecter_proxys(modele, audit_donnees, a)))
proxys <- proxys[order(-proxys$v_cramer), ]
message("Variables du modèle les plus liées aux attributs protégés (V de Cramér) :")
print(head(proxys, 8), row.names = FALSE)
if (any(proxys$v_cramer > 0.3)) message("  ! V de Cramér > 0,3 : justifier la variable ou la retirer.")
message("Lecture : un ratio d'impact < 0,8 n'est pas en soi illégal s'il s'explique par le risque\n",
        "réel (calibration et AUC comparables entre groupes) ; il impose une justification documentée.")


# -----------------------------------------------------------------------------
# 7. Exports
# -----------------------------------------------------------------------------
decision <- cbind(hors_periode[c("id", "date_demande", "montant_demande")],
                  pred$`Hors période`[setdiff(names(pred$`Hors période`), "y")])
decision$decision <- ifelse(decision$score >= seuil, "Accord", "Refus")
decision[decision$decision == "Accord", c("motif_1", "motif_2", "motif_3")] <- NA   # motifs utiles aux refus

write.csv2(modele$grille, file.path(sortie, "grille_score.csv"), row.names = FALSE)
write.csv2(do.call(rbind, lapply(modele$bins, `[[`, "table")), file.path(sortie, "binning_toutes_variables.csv"), row.names = FALSE)
write.csv2(modele$classes, file.path(sortie, "classes_de_risque.csv"), row.names = FALSE)
write.csv2(data.frame(echantillon = rownames(perf), perf), file.path(sortie, "performances.csv"), row.names = FALSE)
write.csv2(bt, file.path(sortie, "backtesting_hors_periode.csv"), row.names = FALSE)
write.csv2(stab, file.path(sortie, "stabilite_psi.csv"), row.names = FALSE)
write.csv2(vc$frequence_selection, file.path(sortie, "validation_croisee_selection.csv"), row.names = FALSE)
write.csv2(challengers, file.path(sortie, "challengers.csv"), row.names = FALSE)
write.csv2(strategie, file.path(sortie, "strategie_seuils.csv"), row.names = FALSE)
write.csv2(equite, file.path(sortie, "audit_equite.csv"), row.names = FALSE)
write.csv2(proxys, file.path(sortie, "audit_proxys.csv"), row.names = FALSE)
write.csv2(decision, file.path(sortie, "decisions_hors_periode.csv"), row.names = FALSE)
graphiques(file.path(sortie, "graphiques.pdf"), modele,
           lapply(pred, function(p) list(pd = p$pd, y = p$y, score = p$score)), vc)
saveRDS(modele, file.path(sortie, "modele_particuliers.rds"))
message("\nRésultats écrits dans ", sortie)

# Utilisation en production :
#   source("fonctions_scoring.R")
#   modele <- readRDS("sorties/particuliers/modele_particuliers.rds")
#   p <- predire_score(modele, nouvelles_demandes)
#   p$decision <- ifelse(p$score >= modele$seuil_acceptation, "Accord", "Refus")
