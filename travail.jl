# ---
# title: Transition végétale
# repository: FelixDeCarufel/Devoir2
# auteurs:
#    - nom: De Carufel
#      prenom: Félix
#      matricule: 20275312
#      github: FelixDeCarufel
#    - nom: Moreau
#      prenom: Maxim
#      matricule: 20269875
#      github: Max80780
#    - nom: Benbezza
#      prenom: Younes
#      matricule: 20315516
#      github: Un-Connu
# ---

# # Introduction

# ## Mise en contexte

# On utilise le concept de succession écologique pour décrire le cycle des espèces présentes sur un territoire suite à une perturbation, qu'elle soit naturelle ou 
# anthropologique. La succession peut être primaire, lorsque la perturbation détruit ou rend indisponible toutes formes de sol, ou secondaire, lorsque la perturbation
# enlève toutes formes de végétaux mais que le sol reste accessible. Lors de succession secondaire, les espèces végétales qui se succèdent chronologiquement sont
# généralement des herbacées, des arbustes puis des arbres tout en ayant des espèces des stades de succession précédent (Ury et al., 2025). 

# Lors de l'aménagement de lignes électriques à haute tension, la présence élevée d'arbres, qui peuvent atteindre de plus grandes hauteurs que les herbacées et arbustes,
# est un enjeu pour la sécurité des infrastructures. L'intervention humaine afin de sélectionner certaines espèces et leur abondance devient alors nécessaire afin 
# qu'elles ne posent pas de problème aux infrastructures lorsque la communauté végétale atteint l'équilibre. De plus, un nombre minimal d'espèces devrait être considéré
# afin de négliger les impacts des modifications anthropologiques sur la biodiversité. Trejo-Pérez et ses collègues (2023) ont prouvé qu'une grande sélection d'espèces
# herbacées avait non seulement un impact positif sur la biodiversité, mais que cela permettait aussi de contrer l'établissement d'arbres plus efficacement.
# ## Question

# Si on devait choisir une espèce d'herbacée et 2 espèces de buissons afin d'aménager un corridor de 200 parcelles sous une ligne à haute tension, lesquelles devraient-on 
# choisir (en se fiant à leur matrice de transition) et à quel ratio devraient-elles être plantées pour que 20% des parcelles soient végétalisées, et que, parmi ces 20%, 
# 30% soient des herbes et 70% soient des buissons tout en s'assurant que la variété de buisson la moins abondante ne représente pas moins de 30% du total des parcelles 
# occupées par des buissons?

# ## Hypothèse et résultats attendus

# L'hypothèse stipule que pour respecter les contraintes, il faudra choisir une matrice de transition qui va favoriser le maintien des parcelles vides, mais qui permet tout de
# même une colonisation modérée par les herbacées et un établissement plus stable et présent des deux espèces de buissons. 
# Comme seulement 20 % des 200 parcelles doivent être végétalisées à l’équilibre, on prévoit que les probabilités de transition vers les états végétalisés 
# demeurent assez faibles, alors que les probabilités de rester dans l’état vide devraient être élevées. On prévoit aussi que l'espèce d'herbacée devra avoir une stabilité plus 
# faible que les buissons pour représenter le 30% des parcelles végétalisées. À l'inverse, les deux espèces de buissons devraient avoir des probabilités, 
# dans la matrice de transition, de persistance plus élevées puisque l’ensemble des buissons doit représenter 70 % de la végétation à l’équilibre. Finalement, les 2 types de 
# buissons devraient avoir des probabilités dans la matrice de transition relativement semblable pour que l’espèce la moins abondante représente au moins 30 % des parcelles 
# occupées par des buissons.

# Pour ce qui est de l'état inital, l'hypothèse stipule qu'il faudra planter probablement un nombre proche du maximum permis de 50 états végétalisés, mais que la majorité de
# ceux-ci devront être des buissons plutôt que des herbacés pour respecter les contraintes. 
# Le modèle déterministe devrait montrer une tendance semblable à celui stochastique, mais sans variabilité aléatoire.

# Alors, les résultats attendus seraient que, dans un modèle de Markov, vu que les probabilités de transition restent constantes dans le temps et que le système est fermé, 
# la succession écologique devrait mener vers une distribution stable des états des parcelles après un certain nombre de générations respectant les contraintes.

# # Description du modèle

# ## États possibles

# Le modèle utilisé est un modèle de succession écologique basé sur une chaine de Markov. Chacune des parcelles du corridor peut passer d'un état à un autre d'une génération à 
# l'autre selon une matrice de transition fixe. Les états possibles sont soit vides (_Barren_) ou végétalisés, l'état végétalisé peut être soit des herbacés (_Grasses_) 
# ou bien deux types de buissons (_Shrubs1_ et _Shrubs2_). 

# ## Corridor et transition d'états

# Le corridor est représenté par 200 parcelles indépendantes qui peuvent se trouver chacune dans un seul état à la fois. Elles peuvent changer d'état à chaque génération selon
# la matrice de transition, qui décrit les probabilités de succession ou bien de persistance des différents états selon l'état de base.

# Le modèle stipule que le système est fermé, donc aucune autre espèce ne peut coloniser le corridor, que les probabilités de transitions sont constantes dans le temps et ne
# dépendent pas de la position des parcelles et du voisinage. De plus, les parcelles sont indépendantes une des autres. 

# ## Simulations

# Il y a deux types de simulation, une stochastique et l'autre déterministe. La stochastique inclut une composante aléatoire dans les transitions entre états et la déterministe
# représente la trajectoire moyenne attendue du système qui elle reste fixe dans le temps. Alors, les probabilités de transition sont constantes et ne dépendent pas de la 
# position spatiale des parcelles ainsi que de la composition du voisinage des parcelles. 
# # Code pour le modèle

# En utilisant autant de sous-sections que nécessaire, expliquez le code que
# vous utilisez pour simuler le modèle. Le texte est aussi important que le code
# en lui-même, et doit faire des liens entre les choix de programmation et la
# question biologique.

# ## Packages nécessaires pour la simulation

using CairoMakie
using Distributions
using ProgressMeter

import Random
Random.seed!(2045)

# ## Conditions initiales

# Vecteur d'états initales des parcelles selon leurs états. On a rajouté un état étant le Shrubs2 pour avoir les 2 types buisson.

s = [150, 0, 25, 25]
states = length(s)
patches = sum(s)

# Matrice de probabilités de transition d'un état à l'autre.

T = zeros(Float64, states, states)
T[1, :] = [85, 6, 7, 7]
T[2, :] = [60, 10, 5, 5]
T[3, :] = [75, 2, 7, 7]
T[4, :] = [80, 6, 7, 7]

# ## Fonction check_tansition_matrix

# Fonction vérifiant que chaque ligne de la matrice de transition correspond à des probabilités. La somme des probabilités sur la ligne de matrice de transition doient 
# être égale à 1 pour que toutes les parcelles soient dans un état quelconque au temps t+1. Si ce n'est pas le cas, la fonction normalise automatiquement les valeurs.

""" 
    function check_tansition_matrix!(T)

Fonction vérifiant que chaque ligne de la matrice de transition correspond à des probabilités. La somme des probabilités sur la ligne de matrice de transition doient 
être égale à 1, si ce n'est pas le cas, la fonction normalise automatiquement les valeurs.

T : La matrice de transition

"""
function check_transition_matrix!(T)
    for ligne in axes(T, 1)
        if sum(T[ligne, :]) != 1
            @warn "La somme de la ligne $(ligne) n'est pas égale à 1 et a été modifiée"
            T[ligne, :] ./= sum(T[ligne, :])
        end
    end
    return T
end

# ## Fonction check_function_arguments

# Cette fonction vérifie que la matrice de transition est bien carrée et que le nombre d’états correspond à la taille de la matrice. Cela permet donc d'éviter des
# incohérences entre les états possibles et leur matrice de transition. 

"""
    function check_function_arguments(transitions, states)

Fonction vérifiant que la matrice de transition est bien carrée et que le nombre d’états correspond à la taille de la matrice. Permet d'éviter des
incohérences entre les états possibles et leur matrice de transition. 

transitions : Matrice des probabilités de transitions
states : Matrice initial de parcelles dans chaque état

"""
function check_function_arguments(transitions, states)
    if size(transitions, 1) != size(transitions, 2)
        throw("La matrice de transition n'est pas carrée")
    end

    if size(transitions, 1) != length(states)
        throw("Le nombre d'états ne correspond pas à la matrice de transition")
    end
    return nothing
end

# ## Fonction verif_nombre_buissons_ini

# On crée une fonction afin de vérfier que le nombre de buissons à l'état initial respecte les conditions imposées.

"""
    verif_nombre_buissons_ini(s)

Vérifier le nombre de buissons à l'état initial (s) et voir s'il respecte les conditions imposées par le devoir.

S'il y a plus de 50 buissons, on les modifie afin qu'ils aient les mêmes proportions, mais avec 50 buissons.
Dans les cas où les nouvelles proportions engendrent une somme de 49, on rajoute une parcelle chez l'état vide.

"""
function verif_nombre_buissons_ini(s)

    ## Si jamais le nombre de buissons dépasse 50, on va:

    if (s[3]+s[4]) > 50

            ## donner un message d'avertissement

            @warn "Il y avait initialement plus que 50 buissons. Les proportions des nombres donnés furent gardées."

            ## et stocker les valeurs qui furent données par l'utilisateur.

            ancienne_valeur3= s[3]
            ancienne_valeur4= s[4]

            ## Par la suite, on change les valeurs, afin qu'elles respectent les conditions, tout en respectant les proportions qui furent données initialement.
            ## Si jamais les nouvelles valeurs donnent des nombres à virgule, on arrondit au nombre le plus bas, puisqu'un buisson et demi n'est pas quelquechose qui est observable dans la réalité.

            s[3]= floor(((ancienne_valeur3 * 50) / (ancienne_valeur3+ancienne_valeur4)))
            s[4]= floor(((ancienne_valeur4 * 50) / (ancienne_valeur3+ancienne_valeur4)))

            ## Dans le cas où les nouvelles proportions, à cause de l'approximation, donnent 49 au lieu de 50, on rajoute une parcelle vide afin qu'il y ait 200 parcelles en tout.

            if (s[3]+s[4]) == 49
                s[1]= s[1]+1
            end

            ## On retourne le nouvel état initial.

            
    end
    return s
end

# ## Fonction verif_etat_initial

# On crée une fonction qui vérifie qu'il n'y a pas d'herbes et qu'il n'y a pas plus de 50 buissons à l'état initial.

"""
    verif_etat_initial(s)

Vérifier l'état initial (s) et voir s'il respecte les conditions imposées par le devoir.

S'il y a des parcelles herbacées, elles sont supprimées.
Si la somme des parcelles n'équivaut pas à 200, les autres valeurs (vides, buisson1 et buisson2) sont modifiées en respectant leurs proportions, mais en étant à 200.

"""
function verif_etat_initial(s)

    ## Si il y a des parcelles herbacées, on rejette cet état initial.

   if s[2] !=0
        @warn"Il ne faut pas qu'il y a d'herbes à l'état initial. Les herbes furent supprimées."
        s[2]=0
    end

    ## Si il y a plus ou moins de 200 parcelles, on rejette aussi cet état initial.

    if sum(s) != 200

        ## On donne un message d'avertissement.

        @warn("Il n'y a pas 200 parcelles.")

        ## On stocke les anciennes valeurs,

        ancienne_valeur1=s[1]
        ancienne_valeur3= s[3]
        ancienne_valeur4= s[4]

        ## On les arrondit au nombre le plus bas, selon des proportions où leur somme est égale à 200.

        s[1]= floor(((ancienne_valeur1 * 200) / (ancienne_valeur3+ancienne_valeur4+ancienne_valeur1)))
        s[3]= floor(((ancienne_valeur3 * 200) / (ancienne_valeur3+ancienne_valeur4+ancienne_valeur1)))
        s[4]= floor(((ancienne_valeur4 * 200) / (ancienne_valeur3+ancienne_valeur4+ancienne_valeur1)))
        
         ## En les arrondissant au nombre le plus bas, il devient possible que la somme des parcelles ne soit pas égale à 200. Si c'est le cas, on rajoute les parcelles manquantes à celles de l'état vide.
        
        if patches == 199
            s[1]=s[1]+1
        end
    end
    return s

end

# Le nombre d'états et le nombre de parcelles dépendent de l'état initial.
# Le nombre d'états équivaut au nombre de conditions possibles à l'état initial, le nombre de parcelles et la somme des nombres de chacun de ses états.

states = length(s)
patches = sum(s)

# ## Fonction _sim_stochastic

# Cette fonction simule le tout de façon stochastique. Pour toutes les parcelles, elle va répartir de façon aléatoire celles-ci vers les états possibles de la génération suivante.
# Tout cela sera fait selon les probabilités de la matrice de transition. Donc, elle représente le caractère aléatoire de la succession écologique simulée ici de façon stochastique.

"""
    function _sim_stochastic!(timeseries, transitions, generation)

Fonction simulant de façon stochastique. Pour toutes les parcelles, elle va répartir de façon aléatoire celles-ci vers les états possibles de la génération suivante.
Tout cela fait selon les probabilités de la matrice de transition.

timeseries : Matrice contenant le nombre de parcelles dans chaque état au fil du temps
transitions : La matrice de transition
generation : Indice de la génération actuelle dans la matrice.
La boucle passe à travers chaque état possible du système

"""
function _sim_stochastic!(timeseries, transitions, generation)
    for state in axes(timeseries, 1)
        pop_change = rand(Multinomial(timeseries[state, generation], transitions[state, :])) 

        ## fait un tirage aléatoire pour savoir comment les parcelles de l’état actuel vont se répartir à la génération suivante.
        ## timeseries[state, generation] = combien de parcelles sont dans cet état en ce moment
        ## transitions[state, :] = les probabilités de passer vers chaque état possible
        ## Multinomial(...) = répartit ces parcelles entre les différents états
        ## rand(...) = fait le tirage au hasard

        timeseries[:, generation+1] .+= pop_change

        ## Cette ligne ajoute le résultat du tirage à la génération suivante.
    
    end
end

# ## Fonction _sim_determ

# Cette fonction simule le tout de façon déterministe. Elle calcule directement la composition attendue des états des différentes parcelles selon 
# l'état de base de celles-ci et la matrice de transition. Donc, à la génération suivante, on va obtenir les états des parcelles en appliquant la matrice de transition. 
# Elle représente donc une tendance moyenne du système, sans effet du hasard.


"""
    function __sim_determ!(timeseries, transitions, generation)

Fonction simulant de façon déterministe. Elle calcule directement la composition attendue des états des différentes parcelles selon l'état de base de celles-ci 
et la matrice de transition. Donc, à la génération suivante, on va obtenir les états des parcelles en appliquant la matrice de transition. Elle représente donc une 
tendance moyenne du système, sans effet du hasard.

timeseries : Matrice contenant le nombre de parcelles dans chaque état au fil du temps
transitions : La matrice de transition
generation : Indice de la génération actuelle dans la matrice.

"""
function _sim_determ!(timeseries, transitions, generation)
    pop_change = (timeseries[:, generation]' * transitions)'

    ## Cette ligne calcule directement combien de parcelles devraient se retrouver dans chaque état à la prochaine génération.
    ## Ici : timeseries[:, generation] = le vecteur des parcelles actuelles
    ## ' = transpose le vecteur pour permettre la multiplication matricielle
    ## * transitions = applique la matrice de transition
    ## le dernier ' remet le résultat en colonne

    timeseries[:, generation+1] .= pop_change
    
    ## place directement les valeurs calculées dans la génération suivante.

end

# ## Fonction simulation

# Cette fonction va exécuter la simulation complète. Elle va initialiser les états des parcelles, puis appliquer les vérifications nécessaires et finalement simuler
# l’évolution du corridor sur plusieurs générations. Elle permet donc d’observer comment la composition végétale va changer au fil du temps jusqu’à l'équilibre.



# transitions → la matrice de transition
# states → le nombre initial de parcelles dans chaque état
# generations → le nombre de générations à simuler (500 par défaut)
# stochastic → permet de choisir une simulation stochastique ou déterministe

"""
    function simulation(transitions, states; generations=500, stochastic=false)

Fonction exécutant la simulation complète. Elle va initialiser les états des parcelles, puis appliquer les vérifications nécessaires et finalement simuler
l'évolution du corridor sur plusieurs générations. Elle permet donc d'observer comment la composition végétale va changer au fil du temps jusqu'à l'équilibre.

transitions : La matrice de transition
states : vecteur contenant le nombre initial de parcelles dans chaque état
generation : Nombre de générations à simuler (500 par défaut)
stochastic : Permet de choisir une simulation stochastique (true) ou déterministe (false)
La fonction retourne la matrice timeseries, qui contient l'évolution du nombre de parcelles dans chaque état au fil du temps.

"""
function simulation(transitions, states; generations=500, stochastic=false)

    _data_type = stochastic ? Int64 : Float32

    ## Cette ligne choisit le type de données utilisé dans la simulation :
    ## Int64 si la simulation est stochastique
    ## Float32 si elle est déterministe

    timeseries = zeros(_data_type, length(states), generations + 1)

    ## crée une matrice qui va enregistrer l’évolution du nombre de parcelles dans chaque état au fil des générations. lignes → les états et colonnes → les générations.

    timeseries[:, 1] = states

    _sim_function! = stochastic ? _sim_stochastic! : _sim_determ!

    ## Cette ligne choisit quelle fonction utiliser pour la simulation :
    ## _sim_stochastic! si on veut une simulation aléatoire
    ## _sim_determ! si on veut une simulation déterministe

    for generation in Base.OneTo(generations) # Répéter pour chaque génération en gros.
        _sim_function!(timeseries, transitions, generation)
    end

    return timeseries # La fonction retourne la matrice timeseries, qui contient l’évolution du nombre de parcelles dans chaque état au fil du temps.
end

# ## Figure générée

# Les noms des états et leurs couleurs dans le graphique.

states_names = ["Barren", "Grasses", "Shrubs1", "Shrubs2"]
states_colors = [:grey40, :orange, :teal, :green]

# Création de la figure et des axes du graphique. L’axe des x représente le nombre de générations et l’axe des y le nombre de parcelles dans chaque état.

f = Figure()
ax = Axis(f[1, 1], xlabel="Nb. générations", ylabel="Nb. parcelles")

# ## Vérification des conditions recherchées lorsqu'on fait des simulations stochastiques et déterministes



"""
    function conditions(transitions, states; gen = 199, iteration = 100, seuil = 0.8)

Fonction qui génère les simulations stochastiques et la simulaton déterministe afin de vérifier si les conditions recherchées dans la question sont respectées. Pour les
simulations stochastiques, un score est enregistré dans la variable condition_sto lorsqu'elles respectent les conditions. Ce score est ensuite comparé à une valeur seuil
afin de savoir si le nombre de simulations stochastiques qui répondent aux critères est assez élevé. Lorsque toutes les conditions sont respectées dans les 2 types de 
simulations, la fonction retourne une figure de l'évolution des parcelles dans le temps ainsi que les états initiales qui ont servis à la générer. Lorsque les conditions 
ne sont pas respectées, le score des simulations stochastiques est fourni.

transitions : La matrice de transition
states : vecteur contenant le nombre initial de parcelles dans chaque état
generation : Nombre de générations à simuler (199 par défaut)
iteration : Le nombre de simulations stochastiques qui sont générées
"""
function conditions(transitions, states; gen = 199, iteration = 200)

    ## S'assurer que les conditions initiales sont respectées

    check_transition_matrix!(transitions)
    check_function_arguments(transitions, states)
    states = verif_nombre_buissons_ini(states)
    states = verif_etat_initial(states)

    patches = sum(states)
    
    ## Vérifier les conditions avec une simulation stochastique

    sto = zeros(Float64, length(states), gen+1, iteration)   # "+1" parce que "_sim_stochastic!" utilise generation+1 pour affecter les valeurs des générations dans timeseries
    condition_sto = 0   # Indicateur du nombre de simulations stochastiques qui respectent les conditions

    ## Réalisation des simulations stochastiques en vérifiant et notant si chaque itération correspond aux critères

    for i in 1:iteration
        sto_sim = simulation(transitions, states; stochastic=true, generations=gen)
        sto[:, :, i] = sto_sim

        for j in eachindex(states)
            lines!(ax, sto_sim[j, :], color=states_colors[j], alpha=0.1)
        end
        
        ## Création de différents objets contenant différents états pour faciliter la compréhension des conditions vérifiées à la fin
        
        final_sto = sto[:, gen+1, i]
        veg_sto = sum(final_sto[2:4])
        shrubs_sto = sum(final_sto[3:4])

        if 0.10 <= veg_sto./patches <= 0.30 && 0.20 <= final_sto[2]./veg_sto <= 0.30 && 0.60 <= shrubs_sto./veg_sto <= 0.80 && min(final_sto[3], final_sto[4])./shrubs_sto >= 0.3
            condition_sto += 1
        end
        
    end

    ## Vérifier les conditions avec une simulation déterministe

    det_sim = simulation(transitions, states; stochastic=false, generations=gen+1)
    for i in eachindex(states)
        lines!(ax, det_sim[i, :], color=states_colors[i], alpha=1, label=states_names[i], linewidth=4)  ## on ne peut pas générer le graph avec les stochastiques seulement, car il faut définir les labels, ce qu'on fait juste avec la déterministe
    end

    ## Création de différents objets contenant différents états pour faciliter la compréhension des conditions vérifiées à la fin
    
    final_det = det_sim[:, end] 
    veg_det = sum(final_det[2:4])
    shrubs_det = sum(final_det[3:4])

    if  0.18 <= veg_det./patches <= 0.22 && 0.28 <= final_det[2]./veg_det <= 0.32 && 0.68 <= shrubs_det./veg_det <= 0.72 && min(final_det[3], final_det[4])./shrubs_det >= 0.3
        condition_det = true
    else
        condition_det = false
    end

    return "Une population initiale de $(states[1]) parcelles vides, $(states[2]) parcelles avec de l'herbe, $(states[3]) parcelles occupées par des buissons de l'espèce 1 et $(states[4]) 
    parcelles occupées par des buissons de l'espèce  avec une matrice de transition de $(transitions), mène à une population finale de $(final_det[1]) parcelles vides, 
    $(final_det[2]) parcelles avec de l'herbe, $(final_det[3]) parcelles avec des buissons de l'espèce 1 et $(final_det[4]) parcelles avec des buissons de l'espèce 2.
    Dans ce scénario, $(condition_sto)% des simulations stochastiques correspondent aux conditions recherchées. Il est $(condition_det) de dire que la simulation 
    déterministe y répond."
    
end

resultat = conditions(T, s)
println(resultat)

# # Résultat des simulations et discussion

# ## Présentation des résultats.

# La figure obtenue à  la fin permet d'observer comment le corridor en dessous de la ligne électrique évolue avec le temps. L'on voit comment les 4 états se 
# stabilisent dans les premières générations afin d'atteindre l'équilibre qui respecte les conditions imposées. L'on voit qu'une majorité des parcelles sont
# vides et que celles-ci sont suivies par les parcelles abritant les 2 espèces des buissons (shrub1 et shrub2), puis par les parcelles herbacées, qui représentent
# la minorité des parcelles. L'équilibre est atteint assez abruptement, toutes les parcelles atteignant se remplissant et devenant herbacées ou buissonées dès la 
# première génération, sans vraiment fluctuer. Ceci est différent dans les simulations stochastiques, qui sont observables avec les lignes pâles derrière les
# lignes principales. Ces simulations montrent d'important fluctuations dû aux éléments aléatoires les régissant. Malgré ces fluctuations, ces simulations suivent
# tout de même l'allure des résultats principaux, leurs proportions étant assez similaires. On voit donc que les résultats de la stimulation déterministe répondent
# aux conditions imposées ce qui est aussi vrai pour la majorité des simulations stochastiques.

axislegend(ax)
tightlimits!(ax)
current_figure()
 
# ## Discussion

# ### Recherche de la matrice de transition idéale

# Lorsque l'on aménage les 200 parcelles au départ, il faut planter 50 buissons et il ne peut pas y avoir d'herbes. L'aménagement
# d'un nombre égal des 2 espèces de buissons augmente les chances d'arriver à des populations qui ont une abondance similaire 
# lorsqu'à l'équilibre. Afin de trouver la meilleure matrice de transition possible pour répondre aux critères demandés, nous avons
# commencé en donnant une matrice de transition identique pour les 4 états qui est de 80% de transition à _barren_ afin de respecter 
# l'obtention de 20% d'espèces végétales, 6% de transition à _grasses_ pour que 30% des parcelles végétalisées soient de l'herbe
# et les 14% de transition restant ont été répartis également entre les 2 espèces de _shrubs_ afin de ne pas avoir une sur-représentation
# d'une des 2 espèces. À partir de cette matrice de départ, des modifications automatiques par essais-erreurs ont été effectuées en gardant en
# tête que les transitions vers l'état _barren_ doivent rester élevé, celles vers _grasses_ doivent être faibles et celles vers les 2
# états de _shrubs_ doivent être similaires.

# ### Limitations du modèle

# Les conditions d'équilibres imposées ont été observées avec un intervalle de plus ou moins 10% dans la simulation de notre modèle stochastique afin de prendre en compte 
# l'effet de la stochasticité. Il est nécessaire d'ajouter un intervalle afin que . Avec celui-ci, 83% des simulations stochastiques répondent aux conditions demandées. Il serait 
# possible d'agrandir cet intervalle afin de pouvoir augmenter encore plus ce pourcentage, mais nous trouvons que plus de 10% commence à diverger grandement des objectifs fixés.

# Les 2 vecteurs de transitions propres aux espèces de buissons, sont quasiment identiques. D'un point de vue biologique, il semble très peu probable d'avoir deux espèces
# différentes qui ont un cycle de succession quasi-identique comme modélisé ici.

# La matrice de transition met une grande emphase sur le retour à l'état _barren_ pour toutes les espèces végétales. Il serait
# possible d'interpréter ces résultats comme une pression anthropologique afin de maintenir l'intégrité des lignes électriques
# ou la présence de nombre élevé d'herbivores. Si la pression des herbivores n'était pas suffisante, introduire une cinquième 
# espèce permettrait de réduire l'entretien requis.

# ### Comparaison au modèle non stochastique

# Dans un modèle non stochastique, le résultat obtenu serait déterminé en ne montrant aucune variation. Ceci ne serait pas trop réaliste, puisque de nombreux processus biologiques
# requièrent de la stochasticité lors de modélisation afin de pouvoir mieux représenter l'élément aléatoire présent dans la réalité. Un modèle déterministe, bien que montrant un 
# résultat "concret" et inchangeable, ne simule pas trop bien les conditions environnementales en constant changement.

# # Références

# Trejo-Pérez, R., Chagnon, P.-L., Boivin, P. et Brisson, J. (2024). Resistance against tree encroachment is driven by richness and identity of herbaceous resident species. Applied Vegetation Science, 27(2), e12781. https://doi.org/10.1111/avsc.12781

# Urry, L. A., Cain, M. L. (Michael L., Wasserman, S. A., Minorsky, P. V., Reece, J. B., Campbell, N. A., Campeau-Péloquin, A., Racine, C. (Enseignante) et Dubé, S. (2020). Biologie (5e édition). Pearson ERPI.