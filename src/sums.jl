type KroneckerSum{T} <: AbstractMatrix{T}
    outer::AbstractMatrix{T}
    inner::AbstractMatrix{T}
end

#ok, evidently there is going to be too much dry violation here.. should correct before proceeding too far.
