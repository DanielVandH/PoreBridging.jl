"""
    get_pore_statistics(mesh)

For the given mesh, returns:

- `pore_area`: The pore area.
- `pore_perimeter`: The pore perimeter.
- `pore_smallest_rad`: The smallest distance from the `centroid` to the pore boundary.
- `centroid`: The centroid of the pore.
"""
function get_pore_statistics(mesh)
    bn = get_boundary_nodes(mesh)
    x_bn = mesh.mesh_information.triangulation.points[1, bn]
    y_bn = mesh.mesh_information.triangulation.points[2, bn]
    push!(x_bn, first(x_bn))
    push!(y_bn, first(y_bn))
    xy_bn = [(x, y) for (x, y) in zip(x_bn, y_bn)]
    xy_mat = [x_bn'; y_bn']
    pore_area = PolygonOps.area(xy_bn)
    pore_perimeter = sum(sqrt((x_bn[i+1] - x_bn[i])^2 + (y_bn[i+1] - y_bn[i])^2) for i in 1:(length(x_bn)-1))
    centroid = PolygonOps.centroid(SVector{2,Float64}.(xy_bn))
    pore_smallest_rad = get_smallest_distance(xy_mat, centroid)
    return pore_area, pore_perimeter, pore_smallest_rad, centroid
end

"""
    get_pore_data(data, saveat_times, type="CRS500"; combine=true)

Returns the pore data for the given `data` set, `saveat_times`, and `type`. If `combine==true`,
then the standard deviation for each summary statistic combines the data from each day in 
`saveat_times`, otherwise it is estimated separately for each day.

Returns:
- `areas`: The void coverages.
- `perimeters`: The normalised void perimeters.
- `radii`: The normalised void radii.
- `σ`: The noise estimates.
"""
function get_pore_data(data, saveat_times, type="CRS500"; combine=true)
    if type == "CRS500" && isfile(normpath(@__DIR__, "..", "data", "square_pore_data.jld2"))
        res = load(joinpath(@__DIR__, "..", "data", "square_pore_data.jld2"))
        areas = res["areas"]
        perimeters = res["perimeters"]
        radii = res["radii"]
        σ = res["σ"]
        return areas, perimeters, radii, σ
    elseif type == "SWV500" && isfile(normpath(@__DIR__, "..", "data", "wave_pore_data.jld2"))
        res = load(joinpath(@__DIR__, "..", "data", "wave_pore_data.jld2"))
        areas = res["areas"]
        perimeters = res["perimeters"]
        radii = res["radii"]
        σ = res["σ"]
        return areas, perimeters, radii, σ
    else # Can't give the original boundary files, sorry. The above imports will load what we did down here.
        ## Define the data 
        geo_types = data[!, :geoType]
        all_crs_idx = findall(geo_types .== type)
        last_crsidx = length(all_crs_idx)
        boundary_list = Vector{Matrix{Float64}}(undef, last_crsidx)
        void_list = Vector{Matrix{Float64}}(undef, last_crsidx)
        day_list = zeros(last_crsidx)
        for (i, j) in pairs(all_crs_idx)
            day = data[j, :Day]
            day_list[i] = day
            image = data[j, :imageIndex]
            pore = data[j, :poreIndex]
            file_name = "data/boundaries/PORE_$(type)DAY$(day)IMAGE$(image)PORE$(pore).dat"
            bnd = readdlm(file_name)
            cx = @views minimum(bnd[:, 1])
            cy = @views minimum(bnd[:, 2])
            @views bnd[:, 1] .-= cx
            @views bnd[:, 2] .-= cy
            boundary_list[i] = bnd
            file_name = "data/boundaries/VOID_$(type)DAY$(day)IMAGE$(image)PORE$(pore).dat"
            bnd = readdlm(file_name)
            @views bnd[:, 1] .-= cx
            @views bnd[:, 2] .-= cy
            void_list[i] = bnd
        end
        ## Now having the boundaries read in, we can compute the statistics
        areas = zeros(last_crsidx)
        perimeters = zeros(last_crsidx)
        radii = zeros(last_crsidx)
        void_areas = zeros(last_crsidx)
        void_perimeters = zeros(last_crsidx)
        void_radii = zeros(last_crsidx)
        for i in eachindex(all_crs_idx)
            bnd = boundary_list[i]
            void = void_list[i]
            A = abs(PolygonOps.area([r for r in eachrow(bnd)]))
            AV = abs(PolygonOps.area([r for r in eachrow(void)]))
            P = sum([norm(bnd[i+1, :] .- bnd[i, :]) for i in 1:(size(bnd, 1)-1)])
            PV = sum([norm(void[i+1, :] .- void[i, :]) for i in 1:(size(void, 1)-1)])
            cxy = PolygonOps.centroid(SVector{2,Float64}.(eachrow(bnd)))
            R = get_smallest_distance(bnd', cxy)
            RV = get_smallest_distance(void', cxy)
            areas[i] = A
            void_areas[i] = AV
            perimeters[i] = P
            void_perimeters[i] = PV
            radii[i] = R
            void_radii[i] = RV
        end
        ## Now scale the statistics
        areas = void_areas ./ areas
        perimeters = void_perimeters ./ perimeters
        radii = void_radii ./ radii
        if length(saveat_times) == 2
            areas = [areas[day_list.==saveat_times[1]], areas[day_list.==saveat_times[2]]]
            perimeters = [perimeters[day_list.==saveat_times[1]], perimeters[day_list.==saveat_times[2]]]
            radii = [radii[day_list.==saveat_times[1]], radii[day_list.==saveat_times[2]]]
        else
            areas = areas[day_list.==saveat_times[1]]
            perimeters = perimeters[day_list.==saveat_times[1]]
            radii = radii[day_list.==saveat_times[1]]
        end
        ## Lastly, we need the noise estimates
        if combine
            sigma_area = std(reduce(vcat, areas)) |> x -> repeat([x], length(saveat_times)) # this just creates [σ_S, σ_S]
            sigma_perimeter = std(reduce(vcat, perimeters)) |> x -> repeat([x], length(saveat_times))
            sigma_radius = std(reduce(vcat, radii)) |> x -> repeat([x], length(saveat_times))
        else
            sigma_area = std.(areas)
            sigma_perimeter = std.(perimeters)
            sigma_radius = std.(radii)
        end
        σ = Dict(:area => sigma_area, :perimeter => sigma_perimeter, :radius => sigma_radius)
        ## Return everything 
        return areas, perimeters, radii, σ
    end
end

"""
    get_data_structs(; save_path=normpath(@__DIR__, "..", "paper", "individual_analysis_figures"), saveat_times=[7, 14])

Returns the `PoreDataOptions` and `ProfileAnalysisSaveOptions` for the profile likelihood analysis.
"""
function get_data_structs(; save_path=normpath(@__DIR__, "..", "paper", "individual_analysis_figures"), saveat_times=[7, 14])
    data = DataFrame(CSV.File(normpath(@__DIR__, "..", "data", "dataset.csv")))
    pore_width = 475.0
    initial_time = 5.0
    final_time = 14.0
    mesh = get_square_pore_geometry(pore_width)
    pore_area, pore_perimeter, pore_smallest_rad, centroid = get_pore_statistics(mesh)
    snapshot_save_times = [5, 7, 14, 25, 28]
    wave_mesh = get_wave_pore_geometry()
    wave_pore_area, wave_pore_perimeter, wave_pore_smallest_rad, wave_centroid = get_pore_statistics(wave_mesh)

    pore_data = PoreDataOptions(;
        data,
        mesh,
        pore_area,
        pore_perimeter,
        pore_smallest_rad,
        centroid,
        wave_mesh,
        wave_pore_area,
        wave_pore_perimeter,
        wave_pore_smallest_rad,
        wave_centroid
    )

    save_options = ProfileAnalysisSaveOptions(;
        save_path,
        initial_time,
        final_time,
        saveat_times,
        snapshot_save_times
    )

    return pore_data, save_options
end