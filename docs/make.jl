using PeriodicInterpolation1D
using Documenter
using Plots

ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(PeriodicInterpolation1D, :DocTestSetup, :(using PeriodicInterpolation1D); recursive = true)

# Add titles of sections and overrides page titles
const titles = Dict(
    # "10-tutorials" => "Tutorials", # example folder title
    "91-developer.md" => "Developer docs",
)

function recursively_list_pages(folder; path_prefix = "")
    pages_list = Any[]
    for file in readdir(folder)
        if file == "index.md"
            # We add index.md separately to make sure it is the first in the list
            continue
        end
        # this is the relative path according to our prefix, not @__DIR__, i.e., relative to `src`
        relpath = joinpath(path_prefix, file)
        # full path of the file
        fullpath = joinpath(folder, relpath)

        if isdir(fullpath)
            # If this is a folder, enter the recursion case
            subsection = recursively_list_pages(fullpath; path_prefix = relpath)

            # Ignore empty folders
            if length(subsection) > 0
                title = if haskey(titles, relpath)
                    titles[relpath]
                else
                    @error "Bad usage: '$relpath' does not have a title set. Fix in 'docs/make.jl'"
                    relpath
                end
                push!(pages_list, title => subsection)
            end

            continue
        end

        if splitext(file)[2] != ".md" # non .md files are ignored
            continue
        elseif haskey(titles, relpath) # case 'title => path'
            push!(pages_list, titles[relpath] => relpath)
        else # case 'title'
            push!(pages_list, relpath)
        end
    end

    return pages_list
end

function list_pages()
    root_dir = joinpath(@__DIR__, "src")
    pages_list = recursively_list_pages(root_dir)

    return ["index.md"; pages_list]
end

makedocs(;
    modules = [PeriodicInterpolation1D],
    authors = "Katharina Kormann <k.kormann@rub.de>, Pierre Navaro <navaro@math.cnrs.fr>, Klaus Reuter <khr@mpcdf.mpg.de> ",
    repo = "https://github.com/juliavlasov/PeriodicInterpolation1D.jl/blob/{commit}{path}#{line}",
    sitename = "PeriodicInterpolation1D.jl",
    format = Documenter.HTML(; canonical = "https://juliavlasov.github.io/PeriodicInterpolation1D.jl"),
    pages = list_pages(),
)

deploydocs(; repo = "github.com/juliavlasov/PeriodicInterpolation1D.jl")
