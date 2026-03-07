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
# qu'elles ne posent pas de problème aux infrastructures lorsque la communauté végétale atteint l'équilibre. De plus, un nombre minimale d'espèce devrait être considérée
# afin de négliger les impacts des modifications anthropologiques sur la biodiversité. Trejo-Pérez et ses collègues (2023) ont prouvés qu'une grande sélection d'espèces
# herbacées avait non seulement un impact positif sur la biodiversité, mais que cela permettait aussi de contrer l'établissement d'arbres plus efficacement.


# ajouter les sources dans la bibliographie : ISBN: 9780135188743 et DOI: 10.1111/avsc.12781

# ## Question

# Si on devait choisir une espèce d'herbacée et 2 espèces de buissons afin d'aménager un corridor de 200 parcelles sous une ligne à haute tension, lesquelles devraient-on 
# choisir (en se fiant à leur matrice de transition) et à quel ratio devraient-elles être plantées pour que 20% des parcelles soient végétalisées, et que, parmi ces 20%, 
# 30% soient des herbes et 70% soient des buissons tout en s'assurant que la variété de buisson la moins abondante ne représente pas moins de 30% du total des parcelles 
# occupées par des buissons?

# ## Hypothèse et résultats attendus

# Clarifiez vos attentes par rapport au résultat de la simulation

# # Description du modèle

# En utilisant autant de sous-sections que nécessaire, expliquez le modèle, ses
# suppositions, et les décisions principales

# Notre modèle suppose que les seuls facteurs qui impactent la composition d'espèces dans notre corridor à la génération "t" est la composition
# en espèce au temps "t-1" et leur matrice de transition. On simule donc un environnement fermé dans lequel aucune autre semence ne peut 
# provenir de l'extérieur. De plus, on assume que le taux d'apparition et de mortalité des espèces présentes est constant de générations en 
# générations. 
# Le modèle utilisé est celui d'un modèle stochastique/déterministe de Markov dans lequel chacune des parcelles peut passer d'un état à un autre selon une matrice de transition
# qui elle reste fixe dans le temps. Alors, les probabilités de transition sont constantes et ne dépendent pas de la position spatiale des parcelles ainsi que de la composition du voisinage des parcelles.

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
# ## Fonction vérifiant que chaque ligne de la matrice de transition correspond à des probabilités.
# La somme des probabilités sur la ligne de matrice de transition doient être égale à 1 pour que toutes les transitions possibles soient bien représentées. Si ce n'est pas le cas,
# la fonction normalise automatiquement les valeurs.

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
        throw("Le nombre d'états ne correspond psa à la matrice de transition")
    end
    return nothing
end

function _sim_stochastic!(timeseries, transitions, generation)
    for state in axes(timeseries, 1)
        pop_change = rand(Multinomial(timeseries[state, generation], transitions[state, :]))
        timeseries[:, generation+1] .+= pop_change
    end
end

function _sim_determ!(timeseries, transitions, generation)
    pop_change = (timeseries[:, generation]' * transitions)'
    timeseries[:, generation+1] .= pop_change
end

function simulation(transitions, states; generations=500, stochastic=false)

    check_transition_matrix!(transitions)
    check_function_arguments(transitions, states)

    _data_type = stochastic ? Int64 : Float32
    timeseries = zeros(_data_type, length(states), generations + 1)
    timeseries[:, 1] = states

    _sim_function! = stochastic ? _sim_stochastic! : _sim_determ!

    for generation in Base.OneTo(generations)
        _sim_function!(timeseries, transitions, generation)
    end

    return timeseries
end

# States
# Barren, Grass, Shrubs1, Shrubs2
s = [100, 0, 0, 0]
states = length(s)
patches = sum(s)

# Transitions
T = zeros(Float64, states, states)
T[1, :] = [110, 8, 0, 2]
T[2, :] = [2, 120, 3, 4]
T[3, :] = [1, 0, 94, 12]
T[4, :] = [12, 0, 4, 15]

states_names = ["Barren", "Grasses", "Shrubs1", "Shrubs2"]
states_colors = [:grey40, :orange, :teal, :blue]

# # Présentation des résultats

# Présentez les résultats des simulations, en faisant un lien avec la question initiale.

f = Figure()
ax = Axis(f[1, 1], xlabel="Nb. générations", ylabel="Nb. parcelles")

# Stochastic simulation
for _ in 1:100
    sto_sim = simulation(T, s; stochastic=true, generations=200)
    for i in eachindex(s)
        lines!(ax, sto_sim[i, :], color=states_colors[i], alpha=0.1)
    end
end

# Deterministic simulation
det_sim = simulation(T, s; stochastic=false, generations=200)
for i in eachindex(s)
    lines!(ax, det_sim[i, :], color=states_colors[i], alpha=1, label=states_names[i], linewidth=4)
end

axislegend(ax)
tightlimits!(ax)
current_figure()
# # Discussion

# Concluez sur le résultat, et sur les limitations du modèle.

# # Comment citer
# On peut aussi citer des références dans le document `references.bib`,
# @ermentrout1993cellular -- la bibliographie sera ajoutée automatiquement à la
# fin du document.
