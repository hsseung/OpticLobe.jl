# combine column assignments and column coordinates to assign cells to p-q coordinates
# `id2pq` - dictionary that assigns cell ID to (p, q) coordinates.
# `id2column` - dictionary that assigns cell ID to column ID.

using DataFrames

# add (p, q) coordinates of cells to table (this will be used to define dictionaries below)
columns_df = hcat(DataFrame(column2pq, [:p, :q]), columns_df)

id2pq = Dict{Int64, Vector{Int64}}()
id2column = Dict{Int64, Int64}()

## dictionaries mapping cell id to pq coordinates, and cell id to column id
for (column, row) in enumerate(eachrow(columns_df))
    for id in skipmissing(row[3:end])
        id2pq[id] = [row.p, row.q]
        id2column[id] = column
    end
end
