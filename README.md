# OpticLobe

Exploring the neurons and connections of the *Drosophila* optic lobe.

[![Build Status](https://github.com/hsseung/OpticLobe.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hsseung/OpticLobe.jl/actions/workflows/CI.yml?query=branch%3Amain)

Code to accompany Matsliah, Yu, et al., 'Neuronal “parts list” and wiring diagram for a visual system'

Version 783 of proofreading.

See https://codex.flywire.ai for more information

## Synopsis

```julia
julia> using OpticLobe

julia> first(sort(Wtt["Tm1", intrinsictypes], rev=true), 10)
10-element Named SparseArrays.SparseVector{Int32, Int32}
celltype  │
──────────┼──────
Pm03      │ 78694
Pm02      │ 70095
Pm08      │ 49431
T3        │ 35806
T2a       │ 35162
Pm05      │ 29965
Tm21      │ 21059
Y3        │ 20896
LMa5      │ 17871
Pm06      │ 17143
```

## Global variables

ind2id
id2ind

intrinsictypes
boundarytypes
centraltypes
visualtypes
alltypes
ind2type 

W
A
Ai

Wtt
Wct
Wtc
infraction
outfraction
inmean
outmean

class2families
family2types

## useful functions

toppre
toppost

showall, strings2ticks, type2ids, convert2arrows

codex_open

ng_open

ng_hyper

## Installation
Building the package downloads data files from the FlyWire Codex, and 
processes the data to create some useful global variables.

