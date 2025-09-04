module Helpers

"""
    bin(df, out_centers; col=1)

Bins a list of events in DataFrame or Matrix `df` into `out_centers` using the specified column (`col`) for the observable.
Returns a vector of counts per bin.
"""
function bin(df, out_edges; col=1)
    counts = zeros(length(out_edges)-1)
    for row in eachrow(df)
        val = row[col]
        # Only bin if val is within the window
        if out_edges[1] <= val < out_edges[end]
            idx = searchsortedlast(out_edges, val)
            if 1 <= idx < length(out_edges)
                counts[idx] += 1
            end
        end
    end
    return counts
end

"""
    rebin(df, out_centers; var_col=1, count_col=2)

Rebins a table of (observable, count) pairs in DataFrame or Matrix `df` into `out_centers`.
`var_col` specifies which column contains the observable to bin in.
`count_col` specifies which column contains the counts.
Returns a vector of summed counts per bin.
"""
function rebin(df, out_edges; var_col=1, count_col=2)
    counts = zeros(length(out_edges)-1)
    for row in eachrow(df)
        val = row[var_col]
        count = row[count_col]
        idx = searchsortedlast(out_edges, val)
        if 1 <= idx < length(out_edges)
            counts[idx] += count
        end
    end
    return counts
end

end # module Helpers