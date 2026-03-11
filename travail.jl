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
#    - nom: Auteur
#      prenom: Troisième
#      matricule: XXXXXXXX
#      github: TroisiAut
# ---

# # Introduction

# Vous pouvez utiliser la syntaxe suivante pour l'_italique_, le **gras**, et
# les caractères `monospace`.

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


# ajouter les sources dans la bibliographie : ISBN: 9780135188743 et DOI: 10.1111/avsc.12781

# ## Question

# Si on devait choisir une espèce d'herbacée et 2 espèces de buissons afin d'aménager un corridor de 200 parcelles sous une ligne à haute tension, lesquelles devraient-on 
# choisir (en se fiant à leur matrice de transition) et à quel ratio devraient-elles être plantées pour que 20% des parcelles soient végétalisées, et que, parmi ces 20%, 
# 30% soient des herbes et 70% soient des buissons tout en s'assurant que la variété de buisson la moins abondante ne représente pas moins de 30% du total des parcelles 
# occupées par des buissons?

# ## Hypothèse et résultats attendus

# Clarifiez vos attentes par rapport au résultat de la simulation
# L'hypothèse stipule que 

# # Description du modèle

# En utilisant autant de sous-sections que nécessaire, expliquez le modèle, ses
# suppositions, et les décisions principales

# Notre modèle suppose que les seuls facteurs qui impactent la composition d'espèces dans notre corridor à la génération "t" est la composition
# en espèce au temps "t-1" et leur matrice de transition. On simule donc un environnement fermé dans lequel aucune autre semence ne peut 
# provenir de l'extérieur. De plus, on assume que le taux d'apparition et de mortalité des espèces présentes est constant de générations en 
# générations. 
# Le modèle utilisé est un modèle stochastique/déterministe de Markov dans lequel chacune des parcelles peut passer d'un état à un autre selon une matrice de transition
# qui elle reste fixe dans le temps. Alors, les probabilités de transition sont constantes et ne dépendent pas de la position spatiale des parcelles ainsi que de la 
# composition du voisinage des parcelles.

# On commence une sous-section avec # ## Titre
# # Code pour le modèle

# En utilisant autant de sous-sections que nécessaire, expliquez le code que
# vous utilisez pour simuler le modèle. Le texte est aussi important que le code
# en lui-même, et doit faire des liens entre les choix de programmation et la
# question biologique.

# ## Packages nécessaires pour la simulation

using CairoMakie
using Distributions

import Random
Random.seed!(2045)
# ## Fonction check_tansition_matrix

# Fonction vérifiant que chaque ligne de la matrice de transition correspond à des probabilités. La somme des probabilités sur la ligne de matrice de transition doient 
# être égale à 1 pour que toutes les parcelles soient dans un état quelconque au temps t+1. Si ce n'est pas le cas, la fonction normalise automatiquement les valeurs.

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

function check_function_arguments(transitions, states)
    if size(transitions, 1) != size(transitions, 2)
        throw("La matrice de transition n'est pas carrée")
    end

    if size(transitions, 1) != length(states)
        throw("Le nombre d'états ne correspond pas à la matrice de transition")
    end
    return nothing
end

# ## Fonction _sim_stochastic

# Cette fonction simule le tout de façon stochastique. Pour toutes les parcelles, elle va répartir de façon aléatoire celles-ci vers les états possibles de la génération suivante.
# Tout cela selon les probabilités de la matrice de transition. Donc, elle représente donc le caractère aléatoire de la succession écologique simuler ici de façon stochastique.

# À partir d'ici les commentaires c'est pour nous on les enlèvera à la fin mais ça nous permet de comprendre le tout:

# timeseries : le tableau qui contient le nombre de parcelles dans chaque état au fil du temps
# transitions : la matrice de transition
# generation : la génération qu’on est en train de simuler
# boucle passe à travers chaque état possible du système

function _sim_stochastic!(timeseries, transitions, generation)
    for state in axes(timeseries, 1)
        pop_change = rand(Multinomial(timeseries[state, generation], transitions[state, :])) 

        # fait un tirage aléatoire pour savoir comment les parcelles de l’état actuel vont se répartir à la génération suivante.
        # timeseries[state, generation] = combien de parcelles sont dans cet état en ce moment
        # transitions[state, :] = les probabilités de passer vers chaque état possible
        # Multinomial(...) = répartit ces parcelles entre les différents états
        # rand(...) = fait le tirage au hasard

        timeseries[:, generation+1] .+= pop_change

        # Cette ligne ajoute le résultat du tirage à la génération suivante.
    
    end
end

# ## Fonction _sim_determ

# Cette fonction simule le tout de façon déterministe. Elle calcule directement la composition attendue des états des différentes parcelles selon 
# l'état de base de celles-ci et la matrice de transition. Donc, à la génération suivante, on va obtenir les états des parcelles en appliquant la matrice de transition. 
# Elle représente donc une tendance moyenne du système, sans effet du hasard.

# À partir d'ici les commentaires c'est pour nous on les enlèvera à la fin mais ça nous permet de comprendre le tout:

function _sim_determ!(timeseries, transitions, generation)
    pop_change = (timeseries[:, generation]' * transitions)'

    # Cette ligne calcule directement combien de parcelles devraient se retrouver dans chaque état à la prochaine génération.
    # Ici : timeseries[:, generation] = le vecteur des parcelles actuelles
    # ' = transpose le vecteur pour permettre la multiplication matricielle
    # * transitions = applique la matrice de transition
    # le dernier ' remet le résultat en colonne

    timeseries[:, generation+1] .= pop_change
    
    # place directement les valeurs calculées dans la génération suivante.

end

# ## Fonction simulation

# Cette fonction va exécuter la simulation complète. Elle va initialiser les états des parcelles, puis appliquer les vérifications nécessaires et finalement simuler
# l’évolution du corridor sur plusieurs générations. Elle permet donc d’observer comment la composition végétale va changer au fil du temps jusqu’à l'équilibre.

# À partir d'ici les commentaires c'est pour nous on les enlèvera à la fin mais ça nous permet de comprendre le tout:

# transitions → la matrice de transition
# states → le nombre initial de parcelles dans chaque état
# generations → le nombre de générations à simuler (500 par défaut)
# stochastic → permet de choisir une simulation stochastique ou déterministe

function simulation(transitions, states; generations=500, stochastic=false)

    check_transition_matrix!(transitions)
    check_function_arguments(transitions, states)

    _data_type = stochastic ? Int64 : Float32

    # Cette ligne choisit le type de données utilisé dans la simulation :
    # Int64 si la simulation est stochastique (on manipule un nombre entier de parcelles)
    # Float32 si elle est déterministe (on peut avoir des valeurs fractionnaires)

    timeseries = zeros(_data_type, length(states), generations + 1)

    # crée une matrice qui va enregistrer l’évolution du nombre de parcelles dans chaque état au fil des générations. lignes → les états et colonnes → les générations.

    timeseries[:, 1] = states

    _sim_function! = stochastic ? _sim_stochastic! : _sim_determ!

    # Cette ligne choisit quelle fonction utiliser pour la simulation :
    # _sim_stochastic! si on veut une simulation aléatoire
    # _sim_determ! si on veut une simulation déterministe

    for generation in Base.OneTo(generations) # Répéter pour chaque génération en gros.
        _sim_function!(timeseries, transitions, generation)
    end

    return timeseries # La fonction retourne la matrice timeseries, qui contient l’évolution du nombre de parcelles dans chaque état au fil du temps.
end

# ## Matrice d'états initales des parcelles selon leurs états

# On a rajouté un état étant le Shrubs2 pour avoir les 2 types buisson.

# States
# Barren, Grass, Shrubs1, Shrubs2
s = [340, 90, 100, 45]
states = length(s)
patches = sum(s)


# On crée une fonction afin de vérfier que le nombre de buissons à l'état initial respecte les conditions imposées.

function verif_nombre_buissons_ini(s)

    # Si jamais le nombre de buissons dépasse 50, on va:

    if (s[3]+s[4]) > 50

            # donner un message d'avertissement

            @warn "Il y avait initialement plus que 50 buissons. Les proportions des nombres donnés furent gardées."

            # et stocker les valeurs qui furent données par l'utilisateur.

            ancienne_valeur3= s[3]
            ancienne_valeur4= s[4]

            # Par la suite, on change les valeurs, afin qu'elles respectent les conditions, tout en respectant les proportions qui furent données initialement.
            # Si jamais les nouvelles valeurs donnent des nombres à virgule, on arrondit au nombre le plus bas, puisqu'un buisson et demi n'est pas quelquechose qui est observable dans la réalité.

            s[3]= floor(((ancienne_valeur3 * 50) / (ancienne_valeur3+ancienne_valeur4)))
            s[4]= floor(((ancienne_valeur4 * 50) / (ancienne_valeur3+ancienne_valeur4)))

            # Dans le cas où les nouvelles proportions, à cause de l'approximation, donnent 49 au lieu de 50, on rajoute une parcelle vide afin qu'il y ait 200 parcelles en tout.

            if (s[3]+s[4]) == 49
                s[1]= s[1]+1
            end

            # On retourne le nouvel état initial.

            return s
    end
end

verif_nombre_buissons_ini(s)

# On crée aussi une fonction qui vérifie l'état initial.

function verif_etat_initial(s)

    # Si il y a des parcelles herbacées, on rejette cet état initial.

   if s[2] !=0
        @warn"Il ne faut pas qu'il y a d'herbes à l'état initial. Les herbes furent supprimées."
        s[2]=0
    end

    # Si il y a plus ou moins de 200 parcelles, on rejette aussi cet état initial.

    if sum(s) != 200
        @warn("Il n'y a pas 200 parcelles.")
        ancienne_valeur1=s[1]
        ancienne_valeur3= s[3]
        ancienne_valeur4= s[4]
        s[1]= floor(((ancienne_valeur1 * 200) / (ancienne_valeur3+ancienne_valeur4+ancienne_valeur1)))
        s[3]= floor(((ancienne_valeur3 * 200) / (ancienne_valeur3+ancienne_valeur4+ancienne_valeur1)))
        s[4]= floor(((ancienne_valeur4 * 200) / (ancienne_valeur3+ancienne_valeur4+ancienne_valeur1)))
        if sum(s)== 199
            s[1]=s[1]+1
        end
    end
    return s
end

verif_etat_initial(s)


# ## Matrice de transition

# Transitions
T = zeros(Float64, states, states)
T[1, :] = [110, 8, 0, 2]
T[2, :] = [2, 120, 3, 4]
T[3, :] = [1, 0, 94, 12]
T[4, :] = [12, 0, 4, 15]

# ## Les noms des états et leurs couleurs dans le graphique.

states_names = ["Barren", "Grasses", "Shrubs1", "Shrubs2"]
states_colors = [:grey40, :orange, :teal, :blue]

# # Présentation des résultats
# ## Figure générée

# Création de la figure et des axes du graphique. L’axe des x représente le nombre de générations et l’axe des y le nombre de parcelles dans chaque état.

f = Figure()
ax = Axis(f[1, 1], xlabel="Nb. générations", ylabel="Nb. parcelles")

 # ## Simulation Stochastique

# Exécutation de la simulation stochastique 100 fois pour voir la variabilité possible de la succession écologique vu le hasard.

for _ in 1:100
    sto_sim = simulation(T, s; stochastic=true, generations=200)
    for i in eachindex(s)
        lines!(ax, sto_sim[i, :], color=states_colors[i], alpha=0.1)
    end
end

# conditions demandées à la fin pour la stochastique
# ## Simulation Deterministe

# Exécution d'une simulation déterministe pour représenter la trajectoire attendue du système quand les probabilités de transition d'états sont 
# appliquées sans variabilité aléatoire. Elle est une ligne noire.

det_sim = simulation(T, s; stochastic=false, generations=200)
for i in eachindex(s)
    lines!(ax, det_sim[i, :], color=states_colors[i], alpha=1, label=states_names[i], linewidth=4)
end

# ## Réorganisation, c'est plus facile pour moi de comprendre comme ça, quand ça va fonctionner je vais mettre le texte et la
# documentation au bon endroit

s = [340, 90, 100, 45]
states = length(s)
patches = sum(s)

T = zeros(Float64, states, states)
T[1, :] = [110, 8, 0, 2]
T[2, :] = [2, 120, 3, 4]
T[3, :] = [1, 0, 94, 12]
T[4, :] = [12, 0, 4, 15]

states_names = ["Barren", "Grasses", "Shrubs1", "Shrubs2"]
states_colors = [:grey40, :orange, :teal, :blue]
f = Figure()
ax = Axis(f[1, 1], xlabel="Nb. générations", ylabel="Nb. parcelles")

function verif_etat_initial(s)

    # Si il y a des parcelles herbacées, on rejette cet état initial.

   if s[2] !=0
        @warn"Il ne faut pas qu'il y a d'herbes à l'état initial. Les herbes furent supprimées."
        s[2]=0
    end

    # Si il y a plus ou moins de 200 parcelles, on rejette aussi cet état initial.

    if sum(s) != 200
        @warn("Il n'y a pas 200 parcelles.")
        ancienne_valeur1=s[1]
        ancienne_valeur3= s[3]
        ancienne_valeur4= s[4]
        s[1]= floor(((ancienne_valeur1 * 200) / (ancienne_valeur3+ancienne_valeur4+ancienne_valeur1)))
        s[3]= floor(((ancienne_valeur3 * 200) / (ancienne_valeur3+ancienne_valeur4+ancienne_valeur1)))
        s[4]= floor(((ancienne_valeur4 * 200) / (ancienne_valeur3+ancienne_valeur4+ancienne_valeur1)))
        if sum(s)== 199
            s[1]=s[1]+1
        end
    end
    return s
end

function verif_nombre_buissons_ini(s)

    # Si jamais le nombre de buissons dépasse 50, on va:

    if (s[3]+s[4]) > 50

            # donner un message d'avertissement

            @warn "Il y avait initialement plus que 50 buissons. Les proportions des nombres donnés furent gardées."

            # et stocker les valeurs qui furent données par l'utilisateur.

            ancienne_valeur3= s[3]
            ancienne_valeur4= s[4]

            # Par la suite, on change les valeurs, afin qu'elles respectent les conditions, tout en respectant les proportions qui furent données initialement.
            # Si jamais les nouvelles valeurs donnent des nombres à virgule, on arrondit au nombre le plus bas, puisqu'un buisson et demi n'est pas quelquechose qui est observable dans la réalité.

            s[3]= floor(((ancienne_valeur3 * 50) / (ancienne_valeur3+ancienne_valeur4)))
            s[4]= floor(((ancienne_valeur4 * 50) / (ancienne_valeur3+ancienne_valeur4)))

            # Dans le cas où les nouvelles proportions, à cause de l'approximation, donnent 49 au lieu de 50, on rajoute une parcelle vide afin qu'il y ait 200 parcelles en tout.

            if (s[3]+s[4]) == 49
                s[1]= s[1]+1
            end

            # On retourne le nouvel état initial.

            return s
    end
end

function check_transition_matrix!(T)
    for ligne in axes(T, 1)
        if sum(T[ligne, :]) != 1
            @warn "La somme de la ligne $(ligne) n'est pas égale à 1 et a été modifiée"
            T[ligne, :] ./= sum(T[ligne, :])
        end
    end
    return T
end

function check_function_arguments(transitions, states)
    if size(transitions, 1) != size(transitions, 2)
        throw("La matrice de transition n'est pas carrée")
    end

    if size(transitions, 1) != length(states)
        throw("Le nombre d'états ne correspond pas à la matrice de transition")
    end
    return nothing
end

function _sim_determ!(timeseries, transitions, generation)
    pop_change = (timeseries[:, generation]' * transitions)'

    # Cette ligne calcule directement combien de parcelles devraient se retrouver dans chaque état à la prochaine génération.
    # Ici : timeseries[:, generation] = le vecteur des parcelles actuelles
    # ' = transpose le vecteur pour permettre la multiplication matricielle
    # * transitions = applique la matrice de transition
    # le dernier ' remet le résultat en colonne

    timeseries[:, generation+1] .= pop_change
    
    # place directement les valeurs calculées dans la génération suivante.

end

function _sim_stochastic!(timeseries, transitions, generation)
    for state in axes(timeseries, 1)
        pop_change = rand(Multinomial(timeseries[state, generation], transitions[state, :])) 

# fait un tirage aléatoire pour savoir comment les parcelles de l’état actuel vont se répartir à la génération suivante.
# timeseries[state, generation] = combien de parcelles sont dans cet état en ce moment
# transitions[state, :] = les probabilités de passer vers chaque état possible
# Multinomial(...) = répartit ces parcelles entre les différents états
# rand(...) = fait le tirage au hasard

        timeseries[:, generation+1] .+= pop_change

# Cette ligne ajoute le résultat du tirage à la génération suivante.
    
    end
end

function simulation(transitions, states; generations=500, stochastic=false)

    check_transition_matrix!(transitions)
    check_function_arguments(transitions, states)
    verif_etat_initial(s)
    verif_nombre_buissons_ini(s)

    _data_type = stochastic ? Int64 : Float32

        # Cette ligne choisit le type de données utilisé dans la simulation :
        # Int64 si la simulation est stochastique (on manipule un nombre entier de parcelles)
        # Float32 si elle est déterministe (on peut avoir des valeurs fractionnaires)

    timeseries = zeros(_data_type, length(states), generations + 1)

        # crée une matrice qui va enregistrer l’évolution du nombre de parcelles dans chaque état au fil des générations. lignes → les états et colonnes → les générations.

    timeseries[:, 1] = states

    _sim_function! = stochastic ? _sim_stochastic! : _sim_determ!

       # Cette ligne choisit quelle fonction utiliser pour la simulation :
       # _sim_stochastic! si on veut une simulation aléatoire
       # _sim_determ! si on veut une simulation déterministe

    for generation in Base.OneTo(generations) # Répéter pour chaque génération en gros.
       _sim_function!(timeseries, transitions, generation)
    end

   return timeseries # La fonction retourne la matrice timeseries, qui contient l’évolution du nombre de parcelles dans chaque état au fil du temps.
end


# J'obtiens "MethodError: no method matching zeros(::Type{Float64}, ::Vector{Int64}, ::Int64, ::Int64)" quand je met "states" dans "sto", mais length(s) qui correspond
# a "states" fonctionne...

function conditions(transitions, states; gen = 199, iteration = 100)
    
    # ## Vérifier les conditions avec une simulation stochastique

    sto = zeros(Float64, length(s), gen+1, iteration)   # "+1" parce que "_sim_stochastic!" utilise generation+1 pour affecter les valeurs des générations dans timeseries
    condition_sto = 0   # Indicateur du nombre de simulations stochastiques qui respectent les conditions

    # Réalisation des simulations stochastiques en vérifiant et notant si chaque itération correspond aux critères

    for i in 1:iteration
        sto_sim = simulation(transitions, states; stochastic=true, generations=gen)
        sto[:, :, i] = sto_sim

        for j in eachindex(states)
            lines!(ax, sto_sim[j, :], color=states_colors[j], alpha=0.1)
        end

        if sum(sto[2:4, 200, i])./patches <= 0.22 && 0.28 <= sto[2, 200, i]./patches <= 0.32 && 0.68 <= sum(sto[3:4, :, i])./patches <= 0.72 && min(sto[3, :, i], sto[4, :, i])./sum(sto[3:4, :, i]) >= 0.3
            condition_sto += 1
        end
        
    end

    # ## Vérifier les conditions avec une simulation déterministe

    det_sim = simulation(T, s; stochastic=false, generations=gen+1)
    for i in eachindex(s)
        lines!(ax, det_sim[i, :], color=states_colors[i], alpha=1, label=states_names[i], linewidth=4)  # on ne peut pas générer le graph savec les stochastiques seulement, car il faut définir les labels, ce qu'on fait juste avec la déterministe
    end
    if  sum(s[2:4])./patches <= 0.22 && 0.28 <= s[2]./patches <= 0.32 && 0.68 <= sum(s[3:4])./patches <= 0.72 && min(s[3], s[4])./sum(s[3:4]) >= 0.3
        condition_det = true
    else
        condition_det = false
    end

    # ## Afficher un graphique si les conditions sont respectées
    if condition_sto/iteration >= 0.8 && condition_det
        axislegend(ax)
        tightlimits!(ax)
        current_figure()
        return(T , s)
    else
        return "$(condition_sto)% des simulations stochastiques correspondent aux conditions recherchées. Il est $(condition_det) de dire que la simulation déterministe y répond"
    end

end
# ## Présentez les résultats des simulations, en faisant un lien avec la question initiale.

# # Discussion

# Concluez sur le résultat, et sur les limitations du modèle.

# # Comment citer
# On peut aussi citer des références dans le document `references.bib`,
# @ermentrout1993cellular -- la bibliographie sera ajoutée automatiquement à la
# fin du document.
