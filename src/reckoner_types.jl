using Distributions

import Tables

abstract type AbstractMatch <: Tables.AbstractRow end

Base.getproperty(m::AbstractMatch, nm::Symbol) = getfield(m, nm)
Base.propertynames(m::AbstractMatch) = fieldnames(typeof(m))
Tables.schema(m::AbstractMatch) = Tables.Schema(propertynames(m), Tuple(typeof.(m)))


abstract type AbstractMatches <: Tables.AbstractColumns end

Base.getproperty(m::AbstractMatches, nm::Symbol) = getfield(m, nm)
Base.propertynames(m::AbstractMatches) = fieldnames(typeof(m))
Tables.schema(m::AbstractMatches) = Tables.Schema(propertynames(m), Tuple(typeof(i[1]) for i in m))