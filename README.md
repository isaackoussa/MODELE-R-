# Modèles de scoring en R : entreprises et particuliers

Deux grilles de score de crédit construites selon la méthode standard des banques (régression logistique sur
*Weight of Evidence*), avec une batterie de contrôles de robustesse. Elles tournent en R de base. `glmnet` et
`ranger` sont optionnels et servent aux modèles challengers.

```bash
Rscript scoring_entreprises.R     # défaillance des entreprises à 1 an
Rscript scoring_particuliers.R    # défaut des emprunteurs particuliers à 12 mois
```

| Fichier | Rôle |
|---------|------|
| `scoring_complet.R` | **tout en un seul fichier** (fonctions + les deux modèles), à copier-coller dans RStudio |
| `tester_une_entreprise.R` | exemple prêt à l'emploi : noter 3 entreprises fictives (saine, moyenne, fragile) |
| `app/` | **application web** : formulaire pour noter une entreprise ou un particulier (ouvrir `app/index.html`) |
| `exporter_modeles_app.R` | met à jour l'application après un réentraînement (écrit `app/modeles.js`) |
| `fonctions_scoring.R` | boîte à outils commune : binning monotone, WoE/IV, sélection, grille de points, classes de risque, validation, PSI, équité |
| `scoring_entreprises.R` | ratios financiers issus du bilan, score, notation de 1 à 8, recalibrage sur un taux de défaut de cycle |
| `scoring_particuliers.R` | score d'octroi, seuil d'acceptation optimisé sur la rentabilité, motifs de refus, audit d'équité |
| `sorties/` | grilles, tables de validation, graphiques (`graphiques.pdf`), modèles sauvegardés (`.rds`) |

## Utiliser vos données

Dans chaque script, renseignez `FICHIER_DONNEES` avec un CSV au format français (séparateur `;`, décimale `,`)
qui contient les colonnes de la fonction `simuler_…()`. Laissé à `NULL`, le script simule un jeu réaliste :
15 000 bilans de 2019 à 2023 pour les entreprises, 20 000 demandes de 2022 à 2024 pour les particuliers. Les
deux jeux contiennent des valeurs manquantes, des valeurs aberrantes et une dérive sur la dernière période.

**Entreprises** : `id, annee, secteur, age, effectif, chiffre_affaires, ebe, total_actif, fonds_propres,
dettes_financieres, tresorerie, resultat_net, actif_circulant, passif_circulant, creances_clients,
ca_precedent, incidents_paiement, defaut`. Les ratios sont calculés par `calculer_ratios()`.

**Particuliers** : `id, date_demande, sexe, nationalite, situation_familiale, age, type_contrat,
revenu_mensuel, anciennete_emploi_mois, statut_logement, nb_credits_en_cours, montant_demande, duree_mois,
anciennete_banque_mois, jours_decouvert_6m, incidents_12m, utilisation_revolving, epargne,
taux_endettement, reste_a_vivre, ratio_montant_revenu, defaut`.

En production :

```r
source("fonctions_scoring.R")
modele <- readRDS("sorties/particuliers/modele_particuliers.rds")
predire_score(modele, nouvelles_demandes)   # pd, score, motif_1..3, classe, pd_classe
```

## Ce qui rend les modèles robustes

| Risque | Parade |
|--------|--------|
| Valeurs extrêmes, erreurs de saisie | les variables sont découpées en classes : une valeur aberrante tombe dans la classe extrême et ne pèse pas plus qu'elle |
| Valeurs manquantes | elles forment une classe à part si elles sont assez nombreuses, sinon elles reçoivent la classe la plus risquée |
| Modalité jamais vue en production | elle reçoit le WoE le plus risqué |
| Sur-apprentissage | chaque classe contient au moins 5 % des observations, avec 7 classes au plus et un taux de défaut monotone imposé |
| Variables redondantes | filtre de corrélation (\|r\| < 0,6) puis VIF < 5 |
| Effets contre-intuitifs | une variable dont le coefficient a le mauvais signe est retirée, de même qu'une variable non significative (p > 5 %) |
| Fuite d'information | alerte quand l'IV dépasse 0,5 |
| Performance optimiste | validation sur un échantillon test et sur une **période postérieure** (*out-of-time*), avec intervalle de confiance bootstrap du Gini |
| Sélection instable | la chaîne entière (binning, sélection, régression) est réapprise sur chacun des 3 × 5 plis de validation croisée, avec la fréquence de sélection de chaque variable |
| Non-linéarités manquées | comparaison avec une forêt aléatoire et un elastic net |
| Dérive de population | PSI du score et de chaque variable |
| PD mal calibrée | Hosmer-Lemeshow, backtesting binomial par classe, alerte automatique et fonction `recalibrer()` |

## Particuliers : conformité et équité

Le scoring de crédit des personnes physiques est classé « à haut risque » par l'AI Act (annexe III), et une
décision entièrement automatisée relève de l'article 22 du RGPD. Le script en tient compte :

- **Attributs protégés exclus du modèle.** `ATTRIBUTS_PROTEGES` contient sexe, nationalité et situation
  familiale. Ces attributs sont conservés uniquement pour l'audit. L'âge reste dans les candidats (il n'est pas
  retenu sur les données simulées) : ajoutez-le à la liste pour l'exclure.
- **Motifs de refus.** Pour chaque dossier, les trois variables qui coûtent le plus de points.
- **Audit d'équité.** Par sexe, nationalité, situation familiale et tranche d'âge : taux d'acceptation, ratio
  d'impact (règle des 4/5), écart de calibration et AUC à l'intérieur de chaque groupe.
- **Détection de proxys.** V de Cramér entre chaque variable du modèle et chaque attribut protégé, avec une alerte
  au-delà de 0,3.

## Résultats sur les données simulées

| | Entreprises | Particuliers |
|---|---|---|
| Variables retenues | capacité de remboursement, autonomie financière, liquidité, trésorerie/actif, incidents de paiement, délai clients, effectif, secteur | jours de découvert, revenu, taux d'endettement, épargne, utilisation du revolving, incidents, contrat, logement |
| Gini test | 55,0 % (IC 95 % : 47–63) | 50,8 % (IC 95 % : 45–56) |
| Gini hors période | 65,0 % | 53,4 % |
| Gini en validation croisée | 58,6 % ± 4,1 | 51,7 % ± 4,0 |
| Challenger forêt aléatoire (test) | 53,6 % | 46,8 % |
| PSI du score | 0,005 (stable) | 0,001 (stable) |

La grille fait au moins aussi bien que la forêt aléatoire : aucune non-linéarité importante ne lui échappe. Sur la
période postérieure, le classement des risques tient, mais le niveau de défaut augmente (choc de 2023 pour les
entreprises, inflation de 2024 pour les particuliers). Le backtesting le détecte et le script recommande de
recalibrer. C'est ce qu'il faut attendre d'un suivi de modèle.

L'audit des particuliers fait ressortir une sous-estimation du risque des 18-25 ans : 11,5 % de défauts observés
contre 6,8 % prédits. Une segmentation ou un recalibrage spécifique est à étudier.

## Application web (`app/`)

Un formulaire qui note une entreprise ou un particulier et affiche le résultat à chaque saisie : probabilité de
défaut, score, classe de risque, ce qui pénalise le plus et détail des points. Le calcul est le même que
`predire_score()` en R, résultat pour résultat.

- **Sur l'ordinateur** : double-cliquer sur `app/index.html`.
- **En ligne** : sur GitHub, *Settings → Pages → Branch : `main`, dossier `/ (root)`*, puis ouvrir
  `https://isaackoussa.github.io/MODELE-R-/app/`.
- **Enregistrer** : donner un nom au dossier puis « Enregistrer ». Les dossiers se rouvrent, se mettent à jour,
  se suppriment et s'exportent en CSV pour Excel. Ouverte depuis l'ordinateur ou GitHub Pages, l'application
  les garde dans le navigateur (sur cet appareil uniquement).
- **Après un réentraînement** sur vos données : `Rscript exporter_modeles_app.R`.
