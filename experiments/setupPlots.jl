function setupErrors(errors)
    correct =filter(err::TestsResults -> abs(err.maxError.err)<0.5, errors)
    faithfull =filter(err::TestsResults -> 1.0>=abs(err.maxError.err)>0.5, errors)
    unfaithfull = filter(err::TestsResults -> abs(err.maxError.err)>1.0, errors)
    return correct, faithfull, unfaithfull
end

function setupBuckets(errors)
    buckets = map(err::TestsResults -> err.buckets, errors)
    return mapreduce(permutedims, hcat, transpose.(buckets))
end

function setupLabels(errors)
    label = map(err::TestsResults ->string(err.f), errors)
    return reshape(label, 1, length(label))
end