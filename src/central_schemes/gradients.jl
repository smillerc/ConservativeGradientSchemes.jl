@inline function curvilinear_cell_gradient(scheme::Central, samples, I::CartesianIndex{N}, centers) where {N}
    return ntuple(Val(N)) do axis
        spacing = logical_spacing(eltype(samples), centers, axis)
        centered_derivative(scheme, samples, I, axis, spacing)
    end
end

@inline function physical_cell_gradient(scheme::Central, samples, I::CartesianIndex{N}, centers, inverse_jacobian) where {N}
    logical = curvilinear_cell_gradient(scheme, samples, I, centers)
    return ntuple(Val(N)) do b
        sum(ntuple(a -> inverse_jacobian[a, b] * logical[a], Val(N)))
    end
end

@inline function curvilinear_cell_gradient!(dest, scheme::Central, samples, I, centers)
    gradient = curvilinear_cell_gradient(scheme, samples, I, centers)
    for d in eachindex(gradient)
        dest[d][I] = gradient[d]
    end
    return nothing
end

@inline function physical_cell_gradient!(dest, scheme::Central, samples, I, centers, inverse_jacobian)
    gradient = physical_cell_gradient(scheme, samples, I, centers, inverse_jacobian)
    for d in eachindex(gradient)
        dest[d][I] = gradient[d]
    end
    return nothing
end
