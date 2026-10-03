# #############################################################################
#  SCORING COMPLET EN R — un seul fichier, à copier-coller dans RStudio
# #############################################################################
#  Partie A : fonctions communes (à exécuter en premier)
#  Partie B : modèle ENTREPRISES  (défaillance à 1 an)
#  Partie C : modèle PARTICULIERS (défaut à 12 mois, audit d'équité)
#
#  Exécution : Ctrl+Maj+Entrée dans RStudio, ou  Rscript scoring_complet.R
#  Résultats : dossier « sorties/ » créé dans le répertoire de travail (getwd()).
#  Packages  : R de base suffit ; install.packages(c("glmnet", "ranger"))
#              active les modèles challengers (facultatif).
# #############################################################################


# #############################################################################
#  PARTIE A — FONCTIONS COMMUNES
# #############################################################################
# =============================================================================
#  fonctions_scoring.R — boîte à outils commune aux deux modèles de scoring
# =============================================================================
#  Méthode : grille de score « à la bâloise »
#    1. découpage en classes (binning) supervisé et MONOTONE de chaque variable ;
#    2. transformation en Weight of Evidence (WoE), valeur d'information (IV) ;
#    3. sélection : IV, corrélations, VIF, signe des coefficients, significativité ;
#    4. régression logistique sur les WoE -> probabilité de défaut (PD) ;
#    5. mise à l'échelle en points (PDO), classes de risque monotones ;
#    6. validation : AUC/Gini, KS, Brier, Hosmer-Lemeshow, bootstrap, validation
#       croisée répétée de TOUTE la chaîne, stabilité (PSI), challenger non linéaire.
#
#  Pourquoi c'est robuste :
#    - le binning rend le modèle insensible aux valeurs extrêmes et traite les
#      valeurs manquantes comme une information (classe propre) ;
#    - la monotonicité imposée empêche le sur-apprentissage des « bosses » ;
#    - toute la chaîne (binning + sélection + régression) est réapprise dans
#      chaque pli de validation croisée : pas de fuite d'information ;
#    - un niveau jamais vu ou une valeur manquante non prévue reçoit la classe
#      la plus risquée (choix prudent) ;
#    - la validation se fait hors échantillon ET hors période (out-of-time).
#
#  Dépendances : R de base (stats, graphics). Optionnelles : glmnet, ranger.
#  Convention : cible = 1 pour un défaut (« mauvais »), 0 sinon (« bon »).
#  WoE = ln(%bons / %mauvais) : un WoE élevé = classe peu risquée.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Découpage de l'échantillon
# -----------------------------------------------------------------------------

#' Découpe stratifiée sur la cible (proportion de défauts identique).
decoupe_stratifiee <- function(y, prop_apprentissage = 0.7) {
  idx <- seq_along(y)
  app <- unlist(lapply(split(idx, y), function(i) i[sample.int(length(i), round(length(i) * prop_apprentissage))]))
  sort(app)
}

#' Plis stratifiés pour la validation croisée.
plis_stratifies <- function(y, k = 5) {
  pli <- integer(length(y))
  for (classe in unique(y)) {
    i <- which(y == classe)
    pli[i] <- sample(rep_len(seq_len(k), length(i)))
  }
  pli
}


# -----------------------------------------------------------------------------
# 2. Binning supervisé monotone + WoE
# -----------------------------------------------------------------------------

.woe <- function(bons, mauvais, B, G) log(((bons + 0.5) / (G + 0.5)) / ((mauvais + 0.5) / (B + 0.5)))
.iv  <- function(bons, mauvais, B, G) {
  pb <- (bons + 0.5) / (G + 0.5); pm <- (mauvais + 0.5) / (B + 0.5)
  (pb - pm) * log(pb / pm)
}

#' Binning d'une variable numérique.
#' Classes fines par quantiles, puis fusions successives jusqu'à ce que :
#'  - chaque classe pèse au moins `part_min` de l'échantillon et contienne
#'    des bons ET des mauvais ;
#'  - le taux de défaut soit monotone (sens donné par la corrélation de Spearman) ;
#'  - il y ait au plus `nb_max` classes.
#' Les manquants forment une classe à part s'ils sont assez nombreux,
#' sinon ils reçoivent le WoE de la classe la plus risquée.
binning_numerique <- function(x, y, nom, nb_initial = 20, part_min = 0.05, nb_max = 8, monotone = TRUE) {
  N <- length(y); B <- sum(y); G <- N - B
  ok <- !is.na(x); xo <- x[ok]; yo <- y[ok]
  valeurs <- sort(unique(xo))
  coupures <- if (length(valeurs) <= nb_initial) valeurs   # variable discrète : une classe par valeur
              else unique(quantile(xo, probs = seq(0, 1, length.out = nb_initial + 1), names = FALSE, type = 7))
  coupures <- coupures[coupures < max(xo)]                 # classes ]-Inf ; c1], ..., ]ck ; +Inf[
  sens <- suppressWarnings(sign(cor(xo, yo, method = "spearman")))
  if (is.na(sens) || sens == 0) sens <- 1

  repeat {
    k <- length(coupures) + 1
    if (k == 1) break
    idx <- findInterval(xo, coupures, left.open = TRUE) + 1
    n <- tabulate(idx, k); m <- tabulate(idx[yo == 1], k); b <- n - m
    tx <- (m + 0.5) / (n + 1)
    # (a) classes trop petites ou pures -> fusion avec le voisin le plus proche
    petites <- which(n < part_min * N | m == 0 | b == 0)
    if (length(petites)) {
      j <- petites[which.min(n[petites])]
      v <- c(j - 1, j + 1); v <- v[v >= 1 & v <= k]
      v <- v[which.min(abs(tx[v] - tx[j]))]
      coupures <- coupures[-min(j, v)]
      next
    }
    # (b) monotonicité -> fusion de la paire en violation la plus proche
    if (monotone) {
      d <- diff(tx)
      viol <- which(d * sens <= 0)
      if (length(viol)) {
        coupures <- coupures[-viol[which.min(abs(d[viol]))]]
        next
      }
    }
    # (c) parcimonie -> au plus nb_max classes (fusion des voisines les plus proches)
    if (k > nb_max) { coupures <- coupures[-which.min(abs(diff(tx)))]; next }
    break
  }

  k <- length(coupures) + 1
  idx <- findInterval(xo, coupures, left.open = TRUE) + 1
  n <- tabulate(idx, k); m <- tabulate(idx[yo == 1], k); b <- n - m
  bornes <- c(-Inf, coupures, Inf)
  libelles <- sprintf("]%s ; %s]", format(signif(bornes[-length(bornes)], 4)), format(signif(bornes[-1], 4)))
  woe <- .woe(b, m, B, G); iv <- .iv(b, m, B, G)

  # valeurs manquantes
  n_na <- sum(!ok); m_na <- sum(y[!ok] == 1)
  if (n_na >= 0.01 * N && m_na > 0 && m_na < n_na) {
    woe_na <- .woe(n_na - m_na, m_na, B, G); iv_na <- .iv(n_na - m_na, m_na, B, G)
    na_propre <- TRUE
  } else {
    woe_na <- min(woe); iv_na <- 0; na_propre <- FALSE   # prudence : classe la plus risquée
  }
  tab <- data.frame(variable = nom, classe = libelles, effectif = n, part = n / N,
                    defauts = m, taux_defaut = m / n, woe = woe, iv = iv)
  if (n_na > 0) tab <- rbind(tab, data.frame(
    variable = nom, classe = if (na_propre) "manquant" else "manquant (-> classe la plus risquée)",
    effectif = n_na, part = n_na / N, defauts = m_na, taux_defaut = m_na / n_na, woe = woe_na, iv = iv_na))
  list(variable = nom, type = "numerique", coupures = coupures, woe = woe, woe_na = woe_na,
       table = tab, iv = sum(tab$iv), sens = sens)
}

#' Binning d'une variable qualitative : modalités triées par taux de défaut,
#' les modalités trop petites sont regroupées avec leur voisine.
#' Une modalité inconnue en production reçoit le WoE le plus risqué.
binning_qualitatif <- function(x, y, nom, part_min = 0.05) {
  N <- length(y); B <- sum(y); G <- N - B
  x <- as.character(x); x[is.na(x)] <- "(manquant)"
  n <- tapply(y, x, length); m <- tapply(y, x, sum)
  ordre <- order((m + 0.5) / (n + 1))
  groupes <- as.list(names(n)[ordre]); n <- n[ordre]; m <- m[ordre]
  repeat {
    k <- length(groupes)
    if (k == 1) break
    b <- n - m
    petites <- which(n < part_min * N | m == 0 | b == 0)
    if (!length(petites)) break
    j <- petites[which.min(n[petites])]
    tx <- (m + 0.5) / (n + 1)
    v <- c(j - 1, j + 1); v <- v[v >= 1 & v <= k]
    v <- v[which.min(abs(tx[v] - tx[j]))]
    a <- min(j, v); z <- max(j, v)
    groupes[[a]] <- c(groupes[[a]], groupes[[z]]); groupes[[z]] <- NULL
    n[a] <- n[a] + n[z]; m[a] <- m[a] + m[z]; n <- n[-z]; m <- m[-z]
  }
  b <- n - m
  woe <- .woe(b, m, B, G); iv <- .iv(b, m, B, G)
  correspondance <- setNames(rep(woe, lengths(groupes)), unlist(groupes))
  tab <- data.frame(variable = nom, classe = vapply(groupes, paste, "", collapse = ", "),
                    effectif = as.vector(n), part = as.vector(n) / N, defauts = as.vector(m),
                    taux_defaut = as.vector(m / n), woe = as.vector(woe), iv = as.vector(iv))
  list(variable = nom, type = "qualitatif", correspondance = correspondance,
       woe_inconnu = min(woe), table = tab, iv = sum(tab$iv))
}

#' Binning de toutes les variables candidates.
binning <- function(donnees, cible, variables, part_min = 0.05, nb_initial = 20, nb_max = 8) {
  y <- donnees[[cible]]
  res <- lapply(variables, function(v) {
    x <- donnees[[v]]
    if (is.numeric(x)) binning_numerique(x, y, v, nb_initial = nb_initial, part_min = part_min, nb_max = nb_max)
    else binning_qualitatif(x, y, v, part_min = part_min)
  })
  setNames(res, variables)
}

#' Transformation WoE d'un jeu de données.
appliquer_woe <- function(donnees, bins) {
  as.data.frame(lapply(bins, function(bin) {
    x <- donnees[[bin$variable]]
    if (bin$type == "numerique") {
      w <- bin$woe[findInterval(x, bin$coupures, left.open = TRUE) + 1]
      w[is.na(x)] <- bin$woe_na
    } else {
      x <- as.character(x); x[is.na(x)] <- "(manquant)"
      w <- unname(bin$correspondance[x])
      w[is.na(w)] <- bin$woe_inconnu
    }
    w
  }), col.names = names(bins))
}


# -----------------------------------------------------------------------------
# 3. Sélection des variables
# -----------------------------------------------------------------------------

vif <- function(X) {
  if (ncol(X) < 2) return(setNames(rep(1, ncol(X)), names(X)))
  sapply(names(X), function(v) {
    r2 <- summary(lm(X[[v]] ~ ., data = X[setdiff(names(X), v)]))$r.squared
    1 / (1 - r2)
  })
}

#' Sélection en entonnoir :
#'  IV >= iv_min  ->  |corrélation| < cor_max (on garde la plus forte IV)
#'  ->  VIF < vif_max  ->  régression : signe attendu (négatif) et p-valeur < p_max.
selectionner_variables <- function(bins, woe, y, iv_min = 0.02, cor_max = 0.6,
                                   vif_max = 5, p_max = 0.05, journal = TRUE) {
  dire <- function(...) if (journal) message(...)
  iv <- sort(sapply(bins, `[[`, "iv"), decreasing = TRUE)
  suspectes <- names(iv)[iv > 0.5]
  if (length(suspectes)) dire("  IV > 0,5 (très prédictives : vérifier l'absence de fuite d'information) : ", paste(suspectes, collapse = ", "))
  retenues <- character(0)
  for (v in names(iv)[iv >= iv_min]) {
    if (!length(retenues) || all(abs(cor(woe[[v]], woe[retenues])) < cor_max)) retenues <- c(retenues, v)
  }
  dire("  IV et corrélations : ", length(retenues), " variables sur ", length(iv))
  repeat {
    vv <- vif(woe[retenues])
    if (max(vv) < vif_max) break
    retenues <- setdiff(retenues, names(which.max(vv)))
  }
  repeat {
    fit <- glm(y ~ ., data = cbind(y = y, woe[retenues]), family = binomial())
    co <- summary(fit)$coefficients[-1, , drop = FALSE]
    mauvais_signe <- rownames(co)[co[, 1] >= 0]
    if (length(mauvais_signe)) {
      retenues <- setdiff(retenues, rownames(co)[which.max(co[, 1])]); next
    }
    if (max(co[, 4]) > p_max && length(retenues) > 1) {
      retenues <- setdiff(retenues, rownames(co)[which.max(co[, 4])]); next
    }
    break
  }
  dire("  Après VIF, signes et significativité : ", length(retenues), " variables")
  list(variables = retenues, iv = iv)
}


# -----------------------------------------------------------------------------
# 4. Ajustement de la chaîne complète et prédiction
# -----------------------------------------------------------------------------

#' Apprend toute la chaîne sur `donnees` et renvoie un objet « modele_score ».
#' parametres : part_min, nb_max, iv_min, cor_max, vif_max, p_max,
#'              score_ref, cote_ref (bons:mauvais au score de référence), pdo.
ajuster_score <- function(donnees, cible, variables, parametres = list(), journal = TRUE) {
  p <- modifyList(list(part_min = 0.05, nb_initial = 20, nb_max = 8, iv_min = 0.02, cor_max = 0.6, vif_max = 5,
                       p_max = 0.05, score_ref = 600, cote_ref = 50, pdo = 20), parametres)
  y <- donnees[[cible]]
  bins <- binning(donnees, cible, variables, part_min = p$part_min, nb_initial = p$nb_initial, nb_max = p$nb_max)
  woe <- appliquer_woe(donnees, bins)
  sel <- selectionner_variables(bins, woe, y, p$iv_min, p$cor_max, p$vif_max, p$p_max, journal)
  fit <- glm(y ~ ., data = cbind(y = y, woe[sel$variables]), family = binomial())

  # Mise à l'échelle : score = offset + facteur * ln(cote bons/mauvais)
  facteur <- p$pdo / log(2)
  offset  <- p$score_ref - facteur * log(p$cote_ref)
  b <- coef(fit); k <- length(sel$variables)
  grille <- do.call(rbind, lapply(sel$variables, function(v) {
    t <- bins[[v]]$table
    t$coefficient <- b[[v]]
    t$points <- round((offset - facteur * b[[1]]) / k - facteur * b[[v]] * t$woe)
    t
  }))
  structure(list(cible = cible, candidates = variables, variables = sel$variables, iv = sel$iv,
                 bins = bins, glm = summary(fit)$coefficients, coefficients = b, facteur = facteur, offset = offset,
                 grille = grille, parametres = p), class = "modele_score")
}

#' Score, PD et motifs (les variables qui coûtent le plus de points) pour de nouvelles données.
predire_score <- function(modele, donnees, nb_motifs = 3) {
  bins <- modele$bins[modele$variables]
  woe <- appliquer_woe(donnees, bins)
  eta <- as.vector(modele$coefficients[1] + as.matrix(woe) %*% modele$coefficients[modele$variables])
  pd <- 1 / (1 + exp(-eta))
  # points par variable (WoE -> points de la grille)
  pts <- sapply(modele$variables, function(v) {
    t <- modele$grille[modele$grille$variable == v, ]
    t$points[match(round(woe[[v]], 10), round(t$woe, 10))]
  })
  pts <- matrix(pts, nrow = nrow(donnees), dimnames = list(NULL, modele$variables))
  max_pts <- sapply(modele$variables, function(v) max(modele$grille$points[modele$grille$variable == v]))
  perte <- sweep(-pts, 2, max_pts, `+`)
  motifs <- t(apply(perte, 1, function(l) {
    o <- order(l, decreasing = TRUE)[seq_len(nb_motifs)]
    ifelse(l[o] > 0, names(l)[o], NA)
  }))
  if (nb_motifs == 1) motifs <- t(motifs)
  res <- data.frame(pd = pd, score = rowSums(pts))
  for (i in seq_len(nb_motifs)) res[[paste0("motif_", i)]] <- motifs[, i]
  if (!is.null(modele$classes)) {
    res$classe <- attribuer_classe(modele$classes, res$score)
    res$pd_classe <- modele$classes$pd[match(res$classe, modele$classes$classe)]
  }
  res
}

#' Recalibrage de la constante pour viser un taux de défaut de long terme
#' (tendance centrale) différent de celui de l'échantillon d'apprentissage.
recalibrer <- function(modele, taux_cible, donnees) {
  eta <- qlogis(predire_score(modele, donnees, 1)$pd)
  delta <- uniroot(function(d) mean(plogis(eta + d)) - taux_cible, c(-10, 10))$root
  modele$coefficients[1] <- modele$coefficients[1] + delta
  k <- length(modele$variables)
  modele$grille$points <- modele$grille$points - round(modele$facteur * delta / k)
  modele
}


# -----------------------------------------------------------------------------
# 5. Classes de risque (échelle de notation)
# -----------------------------------------------------------------------------

#' Découpe le score en classes de taux de défaut strictement croissant
#' (classe 1 = meilleur risque). Départ : quantiles, puis fusion des classes
#' non monotones ou sans défaut.
construire_classes <- function(score, y, nb_classes = 8, effectif_min = 0.02) {
  N <- length(y)
  coupures <- unique(quantile(score, probs = seq(0, 1, length.out = nb_classes + 1), names = FALSE))
  coupures <- coupures[coupures > min(score)]              # classes [c, c'[ : la borne basse est inutile
  repeat {
    k <- length(coupures) + 1
    if (k == 1) break
    i <- findInterval(score, coupures) + 1          # i croissant = meilleur score
    n <- tabulate(i, k); m <- tabulate(i[y == 1], k)
    tx <- (m + 0.5) / (n + 1)
    pb <- which(n < effectif_min * N | m == 0)
    if (length(pb)) {
      j <- pb[1]; v <- if (j == k) j - 1 else j + 1
      coupures <- coupures[-min(j, v)]; next
    }
    d <- diff(tx); viol <- which(d >= 0)            # le taux doit DÉCROÎTRE avec le score
    if (length(viol)) { coupures <- coupures[-viol[which.max(d[viol])]]; next }
    break
  }
  k <- length(coupures) + 1
  i <- findInterval(score, coupures) + 1
  n <- tabulate(i, k); m <- tabulate(i[y == 1], k)
  bornes <- c(-Inf, coupures, Inf)
  # classe 1 = meilleurs scores
  data.frame(classe = k:1, score_min = bornes[-length(bornes)], score_max = bornes[-1],
             effectif = n, defauts = m, pd = m / n)[k:1, ]
}

attribuer_classe <- function(classes, score) {
  cl <- classes[order(classes$score_min), ]
  cl$classe[findInterval(score, cl$score_min[-1]) + 1]
}


# -----------------------------------------------------------------------------
# 6. Mesures de performance et de calibration
# -----------------------------------------------------------------------------

#' AUC (statistique de Mann-Whitney), `risque` croissant avec le défaut.
auc <- function(risque, y) {
  r <- rank(risque); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
gini <- function(risque, y) 2 * auc(risque, y) - 1
ks <- function(risque, y) {
  s <- sort(unique(risque))
  max(abs(ecdf(risque[y == 1])(s) - ecdf(risque[y == 0])(s)))
}
brier <- function(pd, y) mean((pd - y)^2)

#' Test de Hosmer-Lemeshow par déciles de PD.
hosmer_lemeshow <- function(pd, y, g = 10) {
  grp <- cut(pd, unique(quantile(pd, seq(0, 1, length.out = g + 1))), include.lowest = TRUE)
  o <- tapply(y, grp, sum); e <- tapply(pd, grp, sum); n <- tapply(y, grp, length)
  stat <- sum((o - e)^2 / (e * (1 - e / n)))
  ddl <- length(o) - 2
  c(statistique = stat, ddl = ddl, p_valeur = 1 - pchisq(stat, ddl))
}

performances <- function(pd, y) {
  hl <- hosmer_lemeshow(pd, y)
  c(n = length(y), taux_defaut = mean(y), pd_moyenne = mean(pd), auc = auc(pd, y),
    gini = gini(pd, y), ks = ks(pd, y), brier = brier(pd, y), hl_p_valeur = unname(hl["p_valeur"]))
}

#' Signale une PD qui sous-estime le risque observé sur un échantillon
#' (discrimination intacte mais niveau à recaler : voir `recalibrer`).
alerte_calibration <- function(perf, echantillon = "Hors période") {
  x <- perf[echantillon, ]
  if (x[["hl_p_valeur"]] < 0.05 && x[["taux_defaut"]] > x[["pd_moyenne"]])
    message(sprintf(paste0("  -> %s : classement du risque conservé mais niveau sous-estimé ",
                           "(%.1f %% de défauts observés contre %.1f %% prédits) : recalibrer."),
                    echantillon, 100 * x[["taux_defaut"]], 100 * x[["pd_moyenne"]]))
}

#' Intervalle de confiance bootstrap (percentile) du Gini.
gini_bootstrap <- function(pd, y, B = 500, niveau = 0.95) {
  g <- replicate(B, { i <- sample.int(length(y), replace = TRUE); gini(pd[i], y[i]) })
  a <- (1 - niveau) / 2
  c(gini = gini(pd, y), inf = unname(quantile(g, a)), sup = unname(quantile(g, 1 - a)), ecart_type = sd(g))
}

#' Backtesting par classe : taux observé vs PD de la classe, test binomial
#' unilatéral (sous-estimation du risque) et intervalle de Wilson.
backtest_classes <- function(classes, classe_obs, y) {
  do.call(rbind, lapply(classes$classe, function(c) {
    yy <- y[classe_obs == c]; n <- length(yy); d <- sum(yy); pd <- classes$pd[classes$classe == c]
    w <- if (n) prop.test(d, n, correct = FALSE)$conf.int else c(NA, NA)
    data.frame(classe = c, effectif = n, defauts = d, pd_classe = pd,
               taux_observe = if (n) d / n else NA, ic_inf = w[1], ic_sup = w[2],
               p_binomial = if (n) binom.test(d, n, max(pd, 1e-6), alternative = "greater")$p.value else NA)
  }))
}


# -----------------------------------------------------------------------------
# 7. Stabilité
# -----------------------------------------------------------------------------

#' Population Stability Index entre une population de référence et une récente.
#' < 0,10 stable ; 0,10–0,25 à surveiller ; > 0,25 dérive significative.
psi <- function(reference, recent, coupures = NULL, nb = 10) {
  if (is.numeric(reference)) {
    if (is.null(coupures)) coupures <- unique(quantile(reference, seq(0, 1, length.out = nb + 1), na.rm = TRUE, names = FALSE))
    coupures[1] <- -Inf; coupures[length(coupures)] <- Inf
    reference <- cut(reference, coupures, include.lowest = TRUE); recent <- cut(recent, coupures, include.lowest = TRUE)
  }
  niv <- union(levels(factor(reference)), levels(factor(recent)))
  p <- (table(factor(reference, niv)) + 0.5) / (length(reference) + 0.5 * length(niv))
  q <- (table(factor(recent, niv)) + 0.5) / (length(recent) + 0.5 * length(niv))
  sum((q - p) * log(q / p))
}

#' PSI du score et des variables du modèle (sur les WoE, donc sur les classes).
stabilite <- function(modele, reference, recent) {
  s_ref <- predire_score(modele, reference, 1)$score; s_rec <- predire_score(modele, recent, 1)$score
  w_ref <- appliquer_woe(reference, modele$bins[modele$variables])
  w_rec <- appliquer_woe(recent, modele$bins[modele$variables])
  v <- sapply(modele$variables, function(x) psi(factor(w_ref[[x]]), factor(w_rec[[x]])))
  d <- data.frame(element = c("SCORE", modele$variables), psi = c(psi(s_ref, s_rec), v))
  d$diagnostic <- cut(d$psi, c(-Inf, 0.1, 0.25, Inf), labels = c("stable", "à surveiller", "dérive"))
  d
}


# -----------------------------------------------------------------------------
# 8. Validation croisée répétée de toute la chaîne
# -----------------------------------------------------------------------------

#' Réapprend binning + sélection + régression dans chaque pli.
#' Renvoie le Gini par pli et la fréquence de sélection de chaque variable :
#' une variable retenue dans moins de ~60 % des plis est fragile.
validation_croisee <- function(donnees, cible, variables, parametres = list(), k = 5, repetitions = 3) {
  y <- donnees[[cible]]
  res <- list(); selection <- list()
  for (r in seq_len(repetitions)) {
    pli <- plis_stratifies(y, k)
    for (f in seq_len(k)) {
      m <- ajuster_score(donnees[pli != f, ], cible, variables, parametres, journal = FALSE)
      pd <- predire_score(m, donnees[pli == f, ], 1)$pd
      res[[length(res) + 1]] <- data.frame(repetition = r, pli = f,
                                           gini = gini(pd, y[pli == f]), ks = ks(pd, y[pli == f]))
      selection[[length(selection) + 1]] <- m$variables
    }
  }
  res <- do.call(rbind, res)
  freq <- sort(table(factor(unlist(selection), levels = variables)) / length(selection), decreasing = TRUE)
  list(plis = res, gini_moyen = mean(res$gini), gini_ecart_type = sd(res$gini),
       frequence_selection = data.frame(variable = names(freq), frequence = as.vector(freq)))
}


# -----------------------------------------------------------------------------
# 9. Modèles challengers (optionnels)
# -----------------------------------------------------------------------------

imputer <- function(apprentissage, autres) {
  num <- names(apprentissage)[sapply(apprentissage, is.numeric)]
  med <- sapply(apprentissage[num], median, na.rm = TRUE)
  prep <- function(d) {
    for (v in num) { d[[paste0(v, "_na")]] <- as.integer(is.na(d[[v]])); d[[v]][is.na(d[[v]])] <- med[[v]] }
    for (v in setdiff(names(d), num)) if (is.character(d[[v]]) || is.factor(d[[v]])) {
      x <- as.character(d[[v]]); x[is.na(x)] <- "(manquant)"
      d[[v]] <- factor(x, levels = unique(c(sort(unique(as.character(apprentissage[[v]]))), "(manquant)")))
    }
    d
  }
  list(apprentissage = prep(apprentissage), autres = lapply(autres, prep))
}

#' Forêt aléatoire (ranger) sur les variables brutes : si elle fait nettement
#' mieux que la grille, une non-linéarité ou une interaction a été manquée.
challenger_foret <- function(app, autres, cible, variables) {
  if (!requireNamespace("ranger", quietly = TRUE)) { message("  (ranger non installé : challenger ignoré)"); return(NULL) }
  im <- imputer(app[variables], lapply(autres, `[`, variables))
  d <- im$apprentissage; d$.y <- factor(app[[cible]])
  rf <- ranger::ranger(.y ~ ., data = d, probability = TRUE, num.trees = 500, min.node.size = 50,
                       importance = "permutation", seed = 1)
  list(modele = rf, pd = lapply(im$autres, function(x) predict(rf, x)$predictions[, "1"]))
}

#' Régression logistique pénalisée (elastic net) sur les WoE de toutes les
#' variables candidates : contrôle de la stabilité de la sélection.
challenger_glmnet <- function(modele, app, autres) {
  if (!requireNamespace("glmnet", quietly = TRUE)) { message("  (glmnet non installé : challenger ignoré)"); return(NULL) }
  X <- as.matrix(appliquer_woe(app, modele$bins))
  cv <- glmnet::cv.glmnet(X, app[[modele$cible]], family = "binomial", alpha = 0.5, nfolds = 5)
  co <- as.matrix(coef(cv, s = "lambda.1se"))
  list(modele = cv, coefficients = co[co[, 1] != 0, , drop = FALSE],
       pd = lapply(autres, function(d) as.vector(predict(cv, as.matrix(appliquer_woe(d, modele$bins)),
                                                         s = "lambda.1se", type = "response"))))
}


# -----------------------------------------------------------------------------
# 10. Graphiques et rapport
# -----------------------------------------------------------------------------

courbe_roc <- function(pd, y) {
  o <- order(pd, decreasing = TRUE)
  data.frame(fpr = c(0, cumsum(y[o] == 0) / sum(y == 0)), tpr = c(0, cumsum(y[o] == 1) / sum(y == 1)))
}

graphiques <- function(fichier, modele, echantillons, vc = NULL) {
  pdf(fichier, width = 11, height = 8)
  on.exit(dev.off())
  coul <- c("#1f4e79", "#c0504d", "#4f8f3a", "#8064a2")
  # 1. IV
  par(mfrow = c(1, 1), mar = c(5, 14, 3, 2))
  iv <- sort(modele$iv)
  barplot(iv, horiz = TRUE, las = 1, col = ifelse(names(iv) %in% modele$variables, coul[1], "grey80"),
          main = "Valeur d'information (en bleu : variables retenues)", xlab = "IV")
  abline(v = c(0.02, 0.1, 0.3), lty = 2, col = "grey40")
  # 2. ROC
  par(mfrow = c(1, 2), mar = c(5, 5, 3, 2))
  plot(0:1, 0:1, type = "n", xlab = "Taux de faux positifs", ylab = "Taux de vrais positifs", main = "Courbes ROC")
  abline(0, 1, col = "grey")
  for (i in seq_along(echantillons)) {
    r <- courbe_roc(echantillons[[i]]$pd, echantillons[[i]]$y)
    lines(r$fpr, r$tpr, col = coul[i], lwd = 2)
  }
  legend("bottomright", sprintf("%s (Gini %.1f %%)", names(echantillons),
         100 * sapply(echantillons, function(e) gini(e$pd, e$y))), col = coul, lwd = 2, bty = "n")
  # 3. distributions du score
  e <- echantillons[[length(echantillons)]]
  h <- hist(e$score, breaks = 30, plot = FALSE)
  hb <- hist(e$score[e$y == 0], breaks = h$breaks, plot = FALSE); hm <- hist(e$score[e$y == 1], breaks = h$breaks, plot = FALSE)
  plot(hb$mids, hb$density, type = "h", lwd = 6, col = adjustcolor(coul[1], 0.6), ylim = range(0, hb$density, hm$density),
       xlab = "Score", ylab = "Densité", main = paste("Distribution du score -", names(echantillons)[length(echantillons)]))
  lines(hm$mids + diff(h$breaks)[1] / 4, hm$density, type = "h", lwd = 6, col = adjustcolor(coul[2], 0.6))
  legend("topleft", c("sains", "défauts"), col = coul[1:2], lwd = 6, bty = "n")
  # 4. calibration
  par(mfrow = c(1, 2))
  lim <- c(0, 0)
  cal <- lapply(echantillons, function(e) {
    g <- cut(e$pd, unique(quantile(e$pd, seq(0, 1, 0.1))), include.lowest = TRUE)
    data.frame(pd = tapply(e$pd, g, mean), obs = tapply(e$y, g, mean))
  })
  lim <- range(0, unlist(cal))
  plot(lim, lim, type = "n", xlab = "PD prédite moyenne", ylab = "Taux de défaut observé", main = "Calibration par déciles")
  abline(0, 1, col = "grey")
  for (i in seq_along(cal)) points(cal[[i]]$pd, cal[[i]]$obs, col = coul[i], pch = 19, type = "b")
  legend("topleft", names(echantillons), col = coul, pch = 19, bty = "n")
  # 5. classes
  if (!is.null(modele$classes)) {
    bt <- backtest_classes(modele$classes, attribuer_classe(modele$classes, e$score), e$y)
    bp <- barplot(100 * bt$taux_observe, names.arg = bt$classe, col = coul[1], ylim = c(0, 100 * max(bt$ic_sup, na.rm = TRUE)),
                  xlab = "Classe de risque (1 = meilleure)", ylab = "Taux de défaut (%)",
                  main = paste("Défauts par classe -", names(echantillons)[length(echantillons)]))
    arrows(bp, 100 * bt$ic_inf, bp, 100 * bt$ic_sup, angle = 90, code = 3, length = 0.04)
    points(bp, 100 * bt$pd_classe, pch = 18, col = coul[2], cex = 1.6)
    legend("topleft", c("observé (IC 95 %)", "PD de la classe"), pch = c(15, 18), col = coul[1:2], bty = "n")
  }
  # 6. validation croisée
  if (!is.null(vc)) {
    par(mfrow = c(1, 2), mar = c(5, 4, 3, 2))
    boxplot(100 * vc$plis$gini, horizontal = TRUE, col = adjustcolor(coul[1], 0.4),
            main = sprintf("Gini en validation croisée (%d plis)", nrow(vc$plis)), xlab = "Gini (%)")
    fs <- vc$frequence_selection[nrow(vc$frequence_selection):1, ]
    par(mar = c(5, 14, 3, 2))
    barplot(100 * fs$frequence, names.arg = fs$variable, horiz = TRUE, las = 1, col = coul[3],
            main = "Fréquence de sélection (%)", xlim = c(0, 100))
    abline(v = 60, lty = 2)
  }
  invisible(NULL)
}

print.modele_score <- function(x, ...) {
  cat("Grille de score -", length(x$variables), "variables :", paste(x$variables, collapse = ", "), "\n")
  cat(sprintf("Échelle : %d points pour une cote de %d:1, %d points pour doubler la cote\n",
              x$parametres$score_ref, x$parametres$cote_ref, x$parametres$pdo))
  invisible(x)
}


# -----------------------------------------------------------------------------
# 11. Politique d'octroi et équité (scoring des particuliers)
# -----------------------------------------------------------------------------

#' Table de stratégie : pour chaque seuil de score, taux d'acceptation, taux de
#' défaut des acceptés et rentabilité espérée (marge sur les bons, perte en cas
#' de défaut = LGD x exposition sur les mauvais).
table_strategie <- function(score, pd, exposition, marge = 0.06, lgd = 0.45, seuils = NULL, y = NULL) {
  if (is.null(seuils)) seuils <- sort(unique(round(quantile(score, seq(0, 0.6, 0.02), names = FALSE))))
  do.call(rbind, lapply(seuils, function(s) {
    a <- score >= s
    data.frame(seuil = s, taux_acceptation = mean(a),
               pd_moyenne_acceptes = mean(pd[a]),
               taux_defaut_acceptes = if (is.null(y)) NA else mean(y[a]),
               rentabilite_esperee = sum(((1 - pd) * marge - pd * lgd) * exposition * a))
  }))
}

#' Audit d'équité par groupe d'un attribut protégé :
#'  - ratio d'impact = taux d'acceptation du groupe / taux du groupe le plus accepté
#'    (règle des 4/5 : alerte si < 0,8) ;
#'  - écart de calibration (PD moyenne - taux observé) : le score doit être aussi
#'    juste pour chaque groupe ;
#'  - AUC intra-groupe : le score doit classer aussi bien dans chaque groupe.
audit_equite <- function(groupe, accepte, y, pd) {
  groupe <- as.character(groupe); groupe[is.na(groupe)] <- "(manquant)"
  r <- do.call(rbind, lapply(sort(unique(groupe)), function(g) {
    i <- groupe == g
    data.frame(groupe = g, effectif = sum(i), taux_acceptation = mean(accepte[i]),
               taux_defaut_observe = mean(y[i]), pd_moyenne = mean(pd[i]),
               ecart_calibration = mean(pd[i]) - mean(y[i]),
               auc_intra_groupe = if (length(unique(y[i])) == 2) auc(pd[i], y[i]) else NA)
  }))
  r$ratio_impact <- r$taux_acceptation / max(r$taux_acceptation)
  r$alerte <- ifelse(r$ratio_impact < 0.8, "ratio < 0,8", "")
  r
}

#' Détection de variables « proxy » : V de Cramér entre les classes de chaque
#' variable du modèle et un attribut protégé. Au-delà de ~0,3 la variable
#' reconstitue en partie l'attribut et doit être justifiée ou retirée.
detecter_proxys <- function(modele, donnees, attribut) {
  w <- appliquer_woe(donnees, modele$bins[modele$variables])
  a <- as.character(donnees[[attribut]]); a[is.na(a)] <- "(manquant)"
  v <- sapply(modele$variables, function(x) {
    t <- table(w[[x]], a)
    if (min(dim(t)) < 2) return(0)
    chi2 <- suppressWarnings(chisq.test(t, correct = FALSE)$statistic)
    sqrt(unname(chi2) / (sum(t) * (min(dim(t)) - 1)))
  })
  data.frame(attribut = attribut, variable = names(v), v_cramer = round(v, 3))
}


# #############################################################################
#  PARTIE B — MODÈLE ENTREPRISES
# #############################################################################
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
#  Vos données : renseigner FICHIER_DONNEES (CSV ; séparateur « ; », décimale « , »)
#              avec les colonnes décrites dans `simuler_entreprises()`.
#              Sinon un jeu réaliste de 15 000 bilans (2019-2023) est simulé.
# =============================================================================

FICHIER_DONNEES <- NULL          # ex. "mes_bilans.csv"
DOSSIER_SORTIE  <- "sorties/entreprises"
ANNEE_HORS_PERIODE <- 2023       # validation out-of-time : dernière année
TAUX_DEFAUT_LONG_TERME <- NA     # ex. 0.045 : recale les PD sur une moyenne de cycle (NA = pas de recalage)
set.seed(2026)

racine <- getwd()
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
    id = d$id, annee = d$annee, defaut = d$defaut,
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
#   modele <- readRDS("sorties/entreprises/modele_entreprises.rds")
#   predire_score(modele, calculer_ratios(nouveaux_bilans))


# #############################################################################
#  PARTIE C — MODÈLE PARTICULIERS
# #############################################################################
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

racine <- getwd()
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
#   modele <- readRDS("sorties/particuliers/modele_particuliers.rds")
#   p <- predire_score(modele, nouvelles_demandes)
#   p$decision <- ifelse(p$score >= modele$seuil_acceptation, "Accord", "Refus")
