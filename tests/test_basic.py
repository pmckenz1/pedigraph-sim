import recombigraph as rg

def test_simulation_runs():
    gen_list = [
        ["P0", "NA", "NA"],
        ["P1", "NA", "NA"],
        ["F1", "P0", "P1"],
    ]

    model = rg.PedigreeModel(
        pedigree=gen_list,
        chromosomes={"A": 100.0},
        seed=1,
    )

    result = model.simulate()

    assert "F1" in result.individuals
