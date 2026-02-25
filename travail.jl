# ---
# title: Transition végétale
# repository: FelixDeCarufel/Devoir2
# auteurs:
#    - nom: De Carufel
#      prenom: Félix
#      matricule: 20275312
#      github: FelixDeCarufel
#    - nom: Auteur
#      prenom: Deuxième
#      matricule: XXXXXXXX
#      github: DeuxiAut
#    - nom: Auteur
#      prenom: Troisième
#      matricule: XXXXXXXX
#      github: TroisiAut
# ---

# # Introduction

# Le texte de votre introduction va ici.

# Vous pouvez utiliser la syntaxe suivante pour l'_italique_, le **gras**, et
# les caractères `monospace`.

# L'introduction doit avoir les sections suivantes:

# ## Mise en contexte

# Décrivez de manière sommaire la situation biologique

# ## Question

# Introduisez la question que vous allez aborder avec le modèle
# ## Hypothèse et résultats attendus

# Clarifiez vos attentes par rapport au résultat de la simulation

# # Description du modèle

# En utilisant autant de sous-sections que nécessaire, expliquez le modèle, ses
# suppositions, et les décisions principales

# On commence une sous-section avec # ## Titre
# # Code pour le modèle

# En utilisant autant de sous-sections que nécessaire, expliquez le code que
# vous utilisez pour simuler le modèle. Le texte est aussi important que le code
# en lui-même, et doit faire des liens entre les choix de programmation et la
# question biologique.

# ## Packages nécessaires

using CairoMakie
using Distributions

import Random
Random.seed!(2045)
# ## Une autre section

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
# Barren, Grass, Shrubs
s = [100, 0, 0]
states = length(s)
patches = sum(s)

# Transitions
T = zeros(Float64, states, states)
T[1, :] = [110, 8, 0]
T[2, :] = [2, 120, 3]
T[3, :] = [1, 0, 94]
T

states_names = ["Barren", "Grasses", "Shrubs"]
states_colors = [:grey40, :orange, :teal]

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
