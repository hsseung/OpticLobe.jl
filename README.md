# OpticLobe

Exploring the neurons and connections of the *Drosophila* optic lobe.

[![Build Status](https://github.com/hsseung/OpticLobe.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hsseung/OpticLobe.jl/actions/workflows/CI.yml?query=branch%3Amain)

Code to accompany [Matsliah, Yu, et al., "Neuronal parts list and wiring diagram for a visual system"](https://doi.org/10.1038/s41586-024-07981-1) and [Seung, "Predicting visual function by interpreting a neuronal wiring diagram"](https://doi.org/10.1038/s41586-024-07953-5) 

Version 783 of proofreading.

See https://codex.flywire.ai for more information

## Installation

``` julia
julia> import Pkg
julia> Pkg.add(url = "https://github.com/hsseung/OpticLobe.jl")
julia> using OpticLobe
```
You may get "waiting for IO to finish" error messages, which you should be able to ignore.

## Synopsis

IDs of five `Tm1` cells.
``` julia
julia> first(type2ids("Tm1"), 5)
5-element Vector{Int64}:
 720575940599333574
 720575940603884512
 720575940604009062
 720575940604094240
 720575940604151520
```

Top ten target cells of `720575940599333574`, with number of synapses.

``` julia
julia> first(sort(W[Name(720575940599333574), :], rev=true), 10)
10-element Named SparseArrays.SparseVector{Int32, Int32}
cellid             │
───────────────────┼───
720575940620875399 │ 57
720575940629884688 │ 39
720575940616048757 │ 38
720575940615102003 │ 29
720575940611830245 │ 27
720575940622106996 │ 27
720575940619661648 │ 25
720575940629129692 │ 22
720575940639008319 │ 21
720575940626667928 │ 20
```

Cell type of `720575940620875399` (top cell in the list above).
``` julia
julia> ind2type[id2ind(720575940620875399)]
"Pm05"
```

Top ten target cell types of `720575940599333574`, with number of synapses.
``` julia
julia> first(sort(Wct[Name(720575940599333574), :], rev=true), 10)
10-element Named SparseArrays.SparseVector{Int32, Int32}
celltype  │
──────────┼────
Pm03      │ 131
Pm08      │  86
Pm02      │  72
Pm05      │  57
T2a       │  41
LMa5      │  29
T5a       │  29
T3        │  26
Tm21      │  26
Tm4       │  22
```

Number of synapses from cell type `Tm1` to cell type `TmY4`:
```julia
julia> Wtt["Tm1", "TmY4"]
5068
```
Top ten intrinsic cell types targeted by cell type `Tm1`:
``` julia
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

