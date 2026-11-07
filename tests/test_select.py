from structpick import Filters, Ligand, Structure, is_additive, select_structures


def make(pdb_id, chain="A", method="X-ray diffraction", resolution=2.0, coverage=1.0, ligands=()):
    return Structure(pdb_id, chain, method, resolution, coverage, 1, 100, tuple(ligands))


def ids(structures):
    return [s.pdb_id for s in structures]


def test_additives_are_not_bound_ligands():
    s = make(
        "1abc",
        ligands=[
            Ligand("GOL", "glycerol", "A"),
            Ligand("SO4", "", "A"),
            Ligand("HEM", "heme", "A"),
        ],
    )
    assert [lig.code for lig in s.bound_ligands] == ["HEM"]
    assert is_additive("gol")
    assert not is_additive("ATP")


def test_ligands_from_other_chains_are_ignored():
    s = make("1abc", chain="A", ligands=[Ligand("ATP", "", "B")])
    assert s.bound_ligands == ()


def test_resolution_breaks_ties_between_near_full_coverage():
    ranked = select_structures(
        [make("full", resolution=1.5, coverage=1.0), make("near", resolution=1.25, coverage=0.993)]
    )
    assert ids(ranked) == ["near", "full"]


def test_coverage_beats_resolution_when_the_gap_is_real():
    ranked = select_structures(
        [make("partial", resolution=1.0, coverage=0.4), make("full", resolution=2.8)]
    )
    assert ids(ranked) == ["full", "partial"]


def test_entries_without_resolution_go_last():
    ranked = select_structures(
        [make("nmr", method="Solution NMR", resolution=None), make("xray", resolution=3.0)]
    )
    assert ids(ranked) == ["xray", "nmr"]


def test_filters():
    structures = [
        make("good", resolution=1.8, ligands=[Ligand("STI", "imatinib", "A")]),
        make("blurry", resolution=3.5),
        make("partial", coverage=0.3),
        make("cryoem", method="Electron Microscopy", resolution=2.5),
        make("nmr", method="Solution NMR", resolution=None),
        make("apo", resolution=1.9, ligands=[Ligand("GOL", "", "A")]),
    ]
    assert ids(select_structures(structures, Filters(max_resolution=1.9))) == ["good", "apo"]
    assert (
        ids(select_structures(structures, Filters(min_coverage=0.5, method="xray")))[-1] == "blurry"
    )
    assert ids(select_structures(structures, Filters(method="em"))) == ["cryoem"]
    assert ids(select_structures(structures, Filters(require_ligand=True))) == ["good"]
    assert ids(select_structures(structures, Filters(ligand="sti"))) == ["good"]
    assert select_structures(structures, Filters(ligand="ATP")) == []
