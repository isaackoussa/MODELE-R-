# =============================================================================
#  tester_une_entreprise.R — noter vos propres entreprises avec le modèle déjà entraîné
# =============================================================================
#  Pré-requis : avoir lancé une fois scoring_entreprises.R (ou scoring_complet.R),
#  qui crée sorties/entreprises/modele_entreprises.rds.
#  Placez-vous dans le dossier du dépôt : setwd("chemin/vers/MODELE-R-")
#  Montants en milliers d'euros (k€).
# =============================================================================

source("fonctions_scoring.R")
modele <- readRDS("sorties/entreprises/modele_entreprises.rds")

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


# -----------------------------------------------------------------------------
# Les entreprises à tester : modifiez les chiffres ou ajoutez des colonnes
# -----------------------------------------------------------------------------
mes_entreprises <- data.frame(
  id                 = c("BOULANGERIE_SAINE", "GARAGE_MOYEN", "RESTAURANT_FRAGILE"),
  annee              = 2024,
  secteur            = c("Commerce", "Services", "Hôtellerie-restauration"),
  # secteurs possibles : Industrie, Commerce, BTP, Services, Transport, Hôtellerie-restauration
  age                = c(25, 8, 2),          # années d'existence
  effectif           = c(12, 6, 9),          # salariés
  chiffre_affaires   = c(1500, 800, 950),
  ebe                = c(180, 56, -20),      # excédent brut d'exploitation
  total_actif        = c(900, 600, 1100),    # total du bilan
  fonds_propres      = c(450, 150, -40),
  dettes_financieres = c(150, 250, 650),     # emprunts bancaires
  tresorerie         = c(200, 30, 5),
  resultat_net       = c(110, 18, -85),
  actif_circulant    = c(500, 300, 250),     # stocks + créances + trésorerie
  passif_circulant   = c(250, 260, 420),     # dettes fournisseurs, fiscales, sociales
  creances_clients   = c(60, 130, 10),
  ca_precedent       = c(1400, 820, 1050),   # chiffre d'affaires de l'année précédente (NA si inconnu)
  incidents_paiement = c(0, 0, 3)            # incidents de paiement sur 12 mois
)

ratios   <- calculer_ratios(mes_entreprises)
resultat <- cbind(entreprise = ratios$id, predire_score(modele, ratios))
resultat$pd <- sprintf("%.2f %%", 100 * resultat$pd)
print(resultat, row.names = FALSE)

# Lecture :
#   pd        probabilité de faire défaut dans les 12 mois
#   score     points (plus c'est haut, moins c'est risqué ; 600 = 1 défaut pour 50 entreprises)
#   motif_1-3 les ratios qui pénalisent le plus l'entreprise
#   classe    1 = meilleur risque ... 8 = risque le plus élevé
