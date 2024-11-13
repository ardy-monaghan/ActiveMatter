abstract type AbInt end

mutable struct ConcInt <: AbInt
    a::Int64
end

mutable struct DifferentInt{T <: AbInt}
    b::T
    c::Float64
end

a = ConcInt(1)

b = DifferentInt(a, 2.0)