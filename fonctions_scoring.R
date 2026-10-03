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
