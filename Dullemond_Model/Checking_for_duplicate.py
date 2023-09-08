from re import match
from typing import List, TypeVar, Dict

T = TypeVar("T")


l = [ 1,2,3,4,5,6,1,2,5,3,5,6,3,5,6,9

]


def list_in_list(l1: List[T], l2: List[T]) -> bool:
    i = 0
    while i < len(l2) - len(l1):
        if l1 == l2[i : len(l1)]:
            return True
        i += 1
    return False


def test_list_in_list():
    l1 = [1, 2]
    l2 = [1, 2, 3]
    assert list_in_list(l1, l2)
    assert not list_in_list(l2, l1)

    l1 = [2, 3]
    assert list_in_list(l1, l2)
    assert not list_in_list(l2, l2)


def find_duplicating_pattern(list_: List[T], length=3) -> Dict[T, int]:
    matching_patterns_found = {}
    i = 0
    while i + length < len(list_) - 1:
        current_slice = list_[i : i + length]
        haystack = list_[i + length :]
        if list_in_list(current_slice, haystack):
            matching_patterns_found[tuple(current_slice)] = (
                matching_patterns_found.get(tuple(current_slice), 0) + 1
            )
        i += 1
    return matching_patterns_found


def test_find_duplicating_pattern():
    pass


if __name__ == "__main__":
    patterns = find_duplicating_pattern(l)
    print(len(patterns))
    print(patterns)