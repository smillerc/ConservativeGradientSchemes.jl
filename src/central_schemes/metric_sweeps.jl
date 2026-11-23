@kernel function faces_kernel!(dest, scheme, samples, axis, face_domain, low, high)
    i = @index(Global, Linear)
    @inbounds I = face_domain[i]
    @inbounds dest[I] = metric_face(scheme, samples, I, axis, low, high)
end

@kernel function derivative_kernel!(dest, scheme, samples, axis, centers, domain, low, high)
    i = @index(Global, Linear)
    @inbounds I = domain[i]
    low_face = metric_face(scheme, samples, I, axis, low, high)
    high_face = metric_face(scheme, samples, shift(I, axis, 1), axis, low, high)
    spacing = logical_spacing(eltype(samples), centers, axis)
    @inbounds dest[I] = (high_face - low_face) / spacing
end

function stencil_bounds(::Central{Order}, axis, valid_cells) where {Order}
    size(valid_cells, axis) >= Order + 1 ||
        throw(ArgumentError("Central{$Order} requires at least $(Order + 1) valid cells along axis $axis"))
    return first(valid_cells)[axis], last(valid_cells)[axis]
end

function interpolate_metric_edges!(dest, scheme::Central, samples, axis, face_domain::CartesianIndices, valid_cells::CartesianIndices, backend)
    isempty(face_domain) && return nothing
    low, high = stencil_bounds(scheme, axis, valid_cells)
    faces_kernel!(backend)(dest, scheme, samples, axis, face_domain, low, high; ndrange = length(face_domain))
    return nothing
end

function metric_cell_derivative!(dest, scheme::Central, samples, axis, centers, domain::CartesianIndices, valid_cells::CartesianIndices, backend)
    isempty(domain) && return nothing
    low, high = stencil_bounds(scheme, axis, valid_cells)
    derivative_kernel!(backend)(dest, scheme, samples, axis, centers, domain, low, high; ndrange = length(domain))
    return nothing
end

function prepare_first_derivative!(dest, scheme::Central, samples, axis, centers, domain::CartesianIndices, valid_cells::CartesianIndices, backend)
    return metric_cell_derivative!(dest, scheme, samples, axis, centers, domain, valid_cells, backend)
end
