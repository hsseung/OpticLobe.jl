const TARGETDIR = joinpath(@__DIR__, "figures")

# Create target directory if it doesn't exist
if !isdir(TARGETDIR)
    mkdir(TARGETDIR)
end