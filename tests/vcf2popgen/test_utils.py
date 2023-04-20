import pandas as pd

def test_create_pop_dict():
    sample_map = pd.Series(['C', 'A', 'A', 'B'])
    expected_pop_dict = {'C' : 1, 'A': 2, 'B' : 3}
    pop_dict = {}

    for pop_id, pop_label in enumerate(sample_map.unique()):
        pop_dict[pop_label] = pop_id + 1

    assert pop_dict == expected_pop_dict
    