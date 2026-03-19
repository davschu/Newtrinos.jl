@testset "Earth Layers" begin

    @testset "Configuration" begin
        el = Newtrinos.earth_layers.configure()
        @test el isa Newtrinos.earth_layers.EarthLayers
        @test el.cfg isa Newtrinos.earth_layers.PREM
    end

    @testset "ray_circle_path_length" begin
        rcpl = Newtrinos.earth_layers.ray_circle_path_length

        # Ray through center of circle (cz = -1, straight down)
        # For a circle of radius R, detector at y, cz=-1:
        # path length = 2R when detector is at the surface going straight through
        R = 6371.0
        y = R
        cz = -1.0
        L = rcpl(R, y, cz)
        @test L ≈ 2 * R atol=1.0

        # No intersection when radius is too small
        L = rcpl(1.0, 6371.0, -0.1)
        @test L == 0.0

        # Horizontal ray (cz ≈ 0) at the surface
        # Should have short or zero path through small inner spheres
        L_inner = rcpl(3000.0, 6371.0, 0.0)
        @test L_inner == 0.0  # horizontal at surface doesn't reach core

        # Symmetry: path length should be same for ±cz (up vs down geometry)
        # Actually, for upgoing vs downgoing the geometry differs since
        # the detector is at fixed position, but for cz and -cz through full sphere:
        L1 = rcpl(R + 20.0, R - 2.0, -0.5)  # downgoing
        @test L1 > 0.0  # should intersect the atmosphere
    end

    @testset "Layer computation" begin
        el = Newtrinos.earth_layers.configure()
        layers = el.compute_layers()

        # Should have multiple layers
        @test length(layers) > 1

        # Radii should be decreasing (from atmosphere to core)
        @test issorted(layers.radius, rev=true)

        # Densities should be non-negative
        @test all(layers.p_density .>= 0.0)
        @test all(layers.n_density .>= 0.0)

        # Outermost layer should have zero density (atmosphere)
        @test layers.p_density[1] ≈ 0.0
        @test layers.n_density[1] ≈ 0.0

        # Inner layers should have higher density
        @test layers.p_density[end] > layers.p_density[2]
    end

    @testset "Path computation" begin
        el = Newtrinos.earth_layers.configure()
        layers = el.compute_layers()

        # Single coszen value
        paths = el.compute_paths([-0.5], layers)
        @test length(paths) == 1
        @test length(paths[1]) > 0

        # Total path length for straight down should be ~2 * R_earth
        paths_down = el.compute_paths([-1.0], layers)
        total_length = sum(p.length for p in paths_down[1])
        @test total_length ≈ 2 * 6369.0 atol=100.0  # approximately diameter

        # Horizontal path should be shorter than vertical
        paths_horiz = el.compute_paths([-0.1], layers)
        total_horiz = sum(p.length for p in paths_horiz[1])
        @test total_horiz < total_length

        # Multiple coszen values
        czs = [-1.0, -0.5, -0.1]
        paths_multi = el.compute_paths(czs, layers)
        @test length(paths_multi) == 3

        # More vertical = longer path
        lengths = [sum(p.length for p in path) for path in paths_multi]
        @test issorted(lengths, rev=true)
    end
end
