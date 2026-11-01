module NumericalRangeGeometry

export NumericalRange, boundary, crawford_number

using LinearAlgebra: eigvals, tr
using PhaseSectors

include("NumericalRange.jl")
include("utils.jl")

end # module NumericalRangeGeometry
