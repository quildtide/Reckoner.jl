using Distributions

import Tables

abstract type AbstractReckonerMatch <: Tables.AbstractRow end

Base.getproperty(m::AbstractReckonerMatch, nm::Symbol) = getfield(m, nm)
Base.propertynames(m::AbstractReckonerMatch) = fieldnames(typeof(m))
Tables.schema(m::AbstractReckonerMatch) = Tables.Schema(propertynames(m), Tuple(typeof.(m)))


abstract type AbstractReckonerMatches <: Tables.AbstractColumns end

Base.getproperty(m::AbstractReckonerMatches, nm::Symbol) = getfield(m, nm)
Base.propertynames(m::AbstractReckonerMatches) = fieldnames(typeof(m))
Tables.schema(m::AbstractReckonerMatches) = Tables.Schema(propertynames(m), Tuple(typeof(i[1]) for i in m))