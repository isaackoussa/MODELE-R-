# =============================================================================
#  exporter_modeles_app.R — exporte les deux grilles de score pour l'application web
# =============================================================================
#  À relancer après chaque réentraînement : écrit app/modeles.js, lu par app/index.html.
# =============================================================================
library(jsonlite)

exporter <- function(modele) {
  g <- modele$grille
  vars <- lapply(modele$variables, function(v) {
    bin <- modele$bins[[v]]; t <- g[g$variable == v, ]
    pire <- which.min(t$woe)
    if (bin$type == "numerique") {
      k <- length(bin$coupures) + 1
      na_propre <- nrow(t) > k && t$classe[k + 1] == "manquant"
      list(nom = v, type = "num", coef = modele$coefficients[[v]], coupures = I(bin$coupures),
           woe = I(t$woe[1:k]), points = I(t$points[1:k]), libelles = I(trimws(t$classe[1:k])),
           taux = I(t$taux_defaut[1:k]),
           woe_na = if (na_propre) t$woe[k + 1] else t$woe[pire],
           points_na = if (na_propre) t$points[k + 1] else t$points[pire],
           max_points = max(t$points))
    } else {
      corr <- bin$correspondance
      list(nom = v, type = "cat", coef = modele$coefficients[[v]],
           modalites = as.list(setNames(t$points[match(corr, t$woe)], names(corr))),
           woe_modalites = as.list(corr),
           woe_inconnu = bin$woe_inconnu, points_inconnu = t$points[pire],
           max_points = max(t$points))
    }
  })
  list(constante = modele$coefficients[[1]], facteur = modele$facteur, offset = modele$offset,
       variables = vars,
       classes = lapply(seq_len(nrow(modele$classes)), function(i) as.list(modele$classes[i, ])),
       seuil = modele$seuil_acceptation)
}

ent <- readRDS("sorties/entreprises/modele_entreprises.rds")
par <- readRDS("sorties/particuliers/modele_particuliers.rds")
json <- toJSON(list(entreprises = exporter(ent), particuliers = exporter(par)),
               auto_unbox = TRUE, digits = 8, null = "null")
writeLines(c("// Généré par exporter_modeles_app.R — ne pas modifier à la main", paste0("const MODELES = ", json, ";")),
           "app/modeles.js", useBytes = TRUE)
message("app/modeles.js écrit")
