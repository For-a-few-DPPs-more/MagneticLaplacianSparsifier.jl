@testset verbose = true "Magnetic vertex-edge incidence matrix" begin
    g = complete_graph(10)
    @testset "angle=0 recover oriented=$oriented incidence" for oriented in [true]
        graph = MetaGraph(g, :angle, 0.0)
        B = magnetic_incidence_matrix(graph; oriented=oriented)
        B_theo = Graphs.LinAlg.incidence_matrix(graph, eltype(B); oriented=oriented)
        @test B == B_theo
    end
end
