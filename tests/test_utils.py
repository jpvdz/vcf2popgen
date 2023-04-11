import pandas as pd

def test_create_pop_dict(sample_map):
    pop_dict = {}

    for pop_id, pop_label in enumerate(sample_map.unique()):
        pop_dict[pop_label] = pop_id + 1

    print(pop_dict)
    
if __name__ == '__main__':
    test_sample_map = pd.Series(['D', 'A', 'A', 'A', 'B', 'B', 'C', 'C'])
    print(test_sample_map)
    test_create_pop_dict(test_sample_map)
