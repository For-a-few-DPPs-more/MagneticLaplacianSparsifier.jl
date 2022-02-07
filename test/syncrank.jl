
@testset "cumulate angles along tree" begin
    n = 8
    mtsf = MetaGraph(Graph(n))
    t1 = 1
    t2 = 1
    t3 = 1
    t4 = 1
    roots = [1; 6]
    add_edge!(mtsf, 1, 2, :angle, t1)
    add_edge!(mtsf, 2, 3, :angle, t2)
    add_edge!(mtsf, 3, 4, :angle, t3)
    add_edge!(mtsf, 3, 5, :angle, t4)
    add_edge!(mtsf, 6, 7, :angle, t4)
    add_edge!(mtsf, 6, 8, :angle, t4)

    set_prop!(mtsf, :roots, roots)

    angles = cumulate_angles(mtsf)

    result = [0; -1; -2; -3; -3; 0; -1; -1]
    @test norm(angles - result) < 1e-14
end
