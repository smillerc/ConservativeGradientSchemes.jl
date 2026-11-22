@inline sweep_first(::Fluid, q, I, axis, lo, hi) = centered_first(Fluid(), q, I, axis)
@inline sweep_first(tag::Metric, q, I, axis, lo, hi) = directional_first(tag, q, I, axis, lo, hi)
@inline sweep_second(::Fluid, q, first, I, axis, h, lo, hi) =
    interior_curvature(Fluid(), q, first, I, axis, h) / h^2
@inline sweep_second(tag::Metric, q, first, I, axis, h, lo, hi) =
    directional_second(tag, q, first, I, axis, h, lo, hi)

@kernel function first_kernel!(dest, scheme, q, axis, centers, domain, lo, hi)
    i = @index(Global, Linear)
    I = domain[i]
    @inbounds dest[I] = sweep_first(scheme, q, I, axis, lo, hi) / spacing(q[I], centers, axis)
end

@kernel function second_kernel!(dest, scheme, q, first, axis, centers, domain, lo, hi)
    i = @index(Global, Linear)
    I = domain[i]
    @inbounds dest[I] = sweep_second(scheme, q, first, I, axis, spacing(q[I], centers, axis), lo, hi)
end

@inline function metric_face(scheme, q, first, second, I, axis, centers, faces, lo, hi)
    if I[axis] == lo
        return reconstruct_boundary_edge(scheme, q, I, axis, centers, faces, Val(:low), first, second)
    end
    J = shift(I, axis, -1)
    if I[axis] == hi + 1
        return reconstruct_boundary_edge(scheme, q, J, axis, centers, faces, Val(:high), first, second)
    end
    return interior_face(q, J, axis, centers, first, second)
end

@kernel function metric_faces_kernel!(dest, scheme, q, first, second, axis, centers, faces, domain, lo, hi)
    i = @index(Global, Linear)
    I = domain[i]
    @inbounds dest[I] = metric_face(scheme, q, first, second, I, axis, centers, faces, lo, hi)
end

@kernel function metric_derivative_kernel!(dest, scheme, q, first, second, axis, centers, faces, domain, lo, hi)
    i = @index(Global, Linear)
    I = domain[i]
    @inbounds begin
        low = metric_face(scheme, q, first, second, I, axis, centers, faces, lo, hi)
        high = metric_face(scheme, q, first, second, shift(I, axis, 1), axis, centers, faces, lo, hi)
        j = I[axis]
        dest[I] = metric_difference(high, low) / (faces[axis][j + 1] - faces[axis][j])
    end
end

function prepare_first_derivative!(dest, scheme::MEG{6}, q, axis, centers, domain::CartesianIndices, valid_cells::CartesianIndices, backend)
    isempty(domain) && return nothing
    first_kernel!(backend)(dest, Fluid(), q, axis, centers, domain, 0, 0; ndrange = length(domain))
    return nothing
end

function prepare_second_derivative!(dest, scheme::MEG{6}, q, first_derivative, axis, centers, domain::CartesianIndices, valid_cells::CartesianIndices, backend)
    isempty(domain) && return nothing
    second_kernel!(backend)(dest, Fluid(), q, first_derivative, axis, centers, domain, 0, 0; ndrange = length(domain))
    return nothing
end

function prepare_metric_derivatives!(first_derivative, second, scheme::MEG{6}, q, axis, centers, first_domain::CartesianIndices, domain::CartesianIndices, valid_cells::CartesianIndices, backend)
    isempty(first_domain) && isempty(domain) && return nothing
    size(valid_cells, axis) >= 6 || throw(ArgumentError("MEG6 requires at least six valid cells along the derivative axis"))
    lo, hi = first(valid_cells)[axis], last(valid_cells)[axis]
    if !isempty(first_domain)
        first_kernel!(backend)(first_derivative, Metric(), q, axis, centers, first_domain, lo, hi; ndrange = length(first_domain))
    end
    if !isempty(domain)
        second_kernel!(backend)(second, Metric(), q, first_derivative, axis, centers, domain, lo, hi; ndrange = length(domain))
    end
    return nothing
end

function interpolate_metric_edges!(dest, scheme::MEG{6}, q, first_derivative, second, axis, centers, faces, domain::CartesianIndices, valid_cells::CartesianIndices, backend)
    isempty(domain) && return nothing
    metric_faces_kernel!(backend)(dest, scheme, q, first_derivative, second, axis, centers, faces, domain, first(valid_cells)[axis], last(valid_cells)[axis]; ndrange = length(domain))
    return nothing
end

function metric_cell_derivative!(dest, scheme::MEG{6}, q, first_derivative, second, axis, centers, faces, domain::CartesianIndices, valid_cells::CartesianIndices, backend)
    isempty(domain) && return nothing
    metric_derivative_kernel!(backend)(dest, scheme, q, first_derivative, second, axis, centers, faces, domain, first(valid_cells)[axis], last(valid_cells)[axis]; ndrange = length(domain))
    return nothing
end
