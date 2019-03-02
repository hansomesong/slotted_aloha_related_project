# Different from traces in SimuJan20,
# the simulation results are obtained within a 500X500 area.

# % 1 - Slotte avec shadow et macrodiv SC, ex nbPcktOK
# % 2 - Slotte sans shadow et macrodiv SC, ex nbPcktOKSsShadow
# % 3 - Slotte avec shadow sur meilleur, ex nbPcktOKBest
# % 4 - Slotte sans shadow sur meilleur, ex nbPcktOKSsShadowBest
# % 5 - Pur interf moyenne avec shadow et macrodiv SC, ex
# nbPcktOKPureAvgAvecShad
# % 6 - Pur interf moyenne sans shadow et macrodiv SC, ex
# nbPcktOKPureAvgSsShad
# % 7 - Pur interf moyenne avec shadow sur meilleur, ex
# nbPcktOKPureAvgAvecShadBest
# % 8 - Pur interf moyenne  sans shadow sur meilleur, ex
# nbPcktOKPureAvgSsShadBest
# % 9 - Pur interf max avec shadow et macrodiv SC,  ex
# nbPcktOKPureMaxAvecShad
# % 10 - Pur interf max sans shadow sur plus proche, ex nbPcktOKPureMaxSsShad
# % 11 - Pur interf max avec shadow et macrodiv SC, ex
# nbPcktOKPureMaxAvecShadBest
# % 12 - Pur interf max  sans shadow sur meilleur, ex
# nbPcktOKPureMaxSsShadBest
# % 13 - slotte avec shadow et macrodiv MRC,  ex nbPcktMrcOK
# % 14 - Pur interf moyenne avec shadow et macrodiv MRC, ex
# nbPcktMrcOKPureAvg
# % 15 - Pur interf max avec shadow et macrodiv MRC, ex nbPcktMrcOKPureMax

# % 16 - Slotte avec shadow et macrodiv SC mais interferences
# independantes, nbPckScIIDInter
# % 17 - Slotte avec shadow et macrodiv MRC mais interferences
# independantes,nbPckMrcIIDInter
# % 18 - Pure interf moyenne avec shadow et macro SC mais interferences
# independantes
# % 19 - Pure interf moyenne avec shadow et macro MRC mais interferences
# independantes,
# % 20 - Pure interf max avec shadow et macro SC mais interferences
# independantes,
# % 21 - Pure interf max avec shadow et macro MRC mais interferences
# independantes,
# Ce qui est nouveau = à partir de 16 (compris).

# Les numéros 22 à 41 correspondent à MRC limité en slotted Aloha, 42 à 61
# à MRC limité en pure Aloha avg Inteférence, 62 à 81 à MRC limité en pure
# Aloha max Int.

# Chaque ligne correspond à un nombre de branches. Par
# exemple 47 correspond à MRC où on somme les 47-42+1=6 meilleurs SIR en
# pure Aloha avg Interf. Tu remarques que les cas 22, 42 et  62
# correspondent en fait au SC (on ne prend que le meilleur SIR, si cela
# marche, c'est bon ; en revanche si le meilleur est inférieur au seuil,
# tous les autres sont également plus petits).