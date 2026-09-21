#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python version of word-conjugacy.nb (Python 3.8+, standard library only).

A word is a list of integers: i denotes the i-th generator, -i its inverse,
and [] the identity. As in the original program, mg specifies the genus,
and letters range over +/-1, ..., +/-(2*mg).
The reduction rules correspond to the group relation
a1*a2*...*a(2g) = a(2g)*...*a2*a1.

The four original function names are preserved: ListPower, ListPartReplace,
MyReduction, and GetConjugator. Input lists are not modified.
Run this file directly or import its functions into another Python file.

Usage:
    from word_conjugacy_en import GetConjugator, MyReduction
    wd = [1, 2, 3, 4, 6, -1, -2, -4]
    mn, mz = GetConjugator(wd, 3)
    # mn: the selected conjugacy representative; mz: a conjugator.
    # Convention: mn = mz * wd * mz^(-1).

Translation notes:
1. The original replacement order is preserved. After each replacement,
   the scan position advances by exactly one.
2. The original Sort[...][[-1]] selects the last item in sorted order.
   Thus, rotations are maximized in numeric lexicographic order; ties are
   resolved by choosing the largest rotation count. The operation has not
   been changed to min despite the original comment.
3. The default repair=True corrects three issues:
   (a) In MyReduction rules 1-2, the target is Reverse[my]^k2, not my^k2.
       This agrees with the relation a1...a(2g) = a(2g)...a1 used elsewhere.
   (b) An extra list nesting in the final Join in GetConjugator is replaced
       by flat concatenation.
   (c) When the special branch handles tt^q, the reversed candidate also
       retains the q-th power.
4. MyReduction(..., repair=False) preserves the original reduction rules.
   GetConjugator(..., repair=False) preserves both the original reduction
   rules and the original power-branch behavior, but still fixes (b) so the
   result consists of flat integer lists. This mode is for comparison with
   the original notebook; its output should not be relied on to decide
   equality or conjugacy of group elements.
   The original program compares lists by length before comparing their
   entries; compatibility mode preserves this behavior.
5. A concrete example of the original rule issue (mg=2):
   MyReduction([1,2,3,4,-1,-2,-3,-4], 2, repair=False)
   returns [1,3,4,-1,-3,-4]. The input is the relator for the relation above
   and should represent the identity. The default repaired version returns [].
6. MyReduction and GetConjugator still follow the notebook's overall
   algorithm. Translating and testing the program does not constitute a
   mathematical proof that it produces a unique conjugacy normal form
   for every input.

Reference for sorting semantics:
https://reference.wolfram.com/language/ref/Sort.html
"""

from numbers import Integral
from typing import Iterable, List, Optional, Sequence, Tuple


Word = List[int]

__all__ = [
    "ListPower", "ListPartReplace", "MyReduction", "GetConjugator",
    "list_power", "list_part_replace", "my_reduction", "get_conjugator",
]


def _integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError("{} must be an integer".format(name))
    return int(value)


def _integers(values: Iterable[int], name: str) -> Word:
    return [_integer(x, name + " element") for x in values]


def _checked_word(values: Iterable[int], mg: int) -> Tuple[Word, int]:
    mg = _integer(mg, "mg")
    if mg < 1:
        raise ValueError("mg must be a positive integer")
    word = _integers(values, "word")
    if any(x == 0 or abs(x) > 2 * mg for x in word):
        raise ValueError("Letters must belong to +/-1, ..., +/-{}".format(2 * mg))
    return word, mg


def _inverse(word: Sequence[int]) -> Word:
    """Invert a word by reversing its order and negating each letter."""
    return [-x for x in reversed(word)]


def _rotate_left(word: Sequence[int], count: int) -> Word:
    if not word:
        return []
    count %= len(word)
    return list(word[count:]) + list(word[:count])


def _common_prefix_length(left: Sequence[int], right: Sequence[int]) -> int:
    count = 0
    for a, b in zip(left, right):
        if a != b:
            break
        count += 1
    return count


def _max_rotation(word: Sequence[int]) -> Tuple[Word, int]:
    """Implement the original Sort[Table[Append[RotateLeft[...], {i}], ...]][[-1]]."""
    if not word:
        raise ValueError("The original algorithm did not extract a nonempty core word; no rotation can be selected")
    # Rotation counts are 1,...,len(word), not 0,...,len(word)-1.
    # The second tuple item reproduces the original tie-breaking rule.
    return max((_rotate_left(word, i), i) for i in range(1, len(word) + 1))


def ListPower(lt: Iterable[int], pr: int) -> Word:
    """Concatenate copies as in ListPower[lt, pr]; invert for pr < 0; return [] for pr = 0."""
    lt = _integers(lt, "lt")
    pr = _integer(pr, "pr")
    nl = lt * abs(pr)
    return _inverse(nl) if pr < 0 else nl


def ListPartReplace(li: Iterable[int], fr: Iterable[int], to: Iterable[int]) -> Word:
    """Scan left to right and replace contiguous subwords, as in ListPartReplace[li, fr, to].

    As in the original For loop, advance the index by 1 after a replacement
    and use the updated list length to determine when to stop.
    This function performs one scan; MyReduction repeats the reduction passes.
    fr must be nonempty. The reduction rules never increase word length.
    """
    st, fr, to = list(li), list(fr), list(to)
    if not fr:
        raise ValueError("fr must not be an empty list")
    k1 = 0
    while k1 <= len(st) - len(fr):
        if st[k1:k1 + len(fr)] == fr:
            st[k1:k1 + len(fr)] = to
        k1 += 1
    return st


def MyReduction(li: Iterable[int], mg: int, *, repair: bool = True) -> Word:
    """Reduce repeatedly in the notebook's rule order, as in MyReduction[li, mg].

    repair=True corrects the target word in rules 1-2; False preserves
    the original rules for direct comparison.
    """
    st, mg = _checked_word(li, mg)
    if not isinstance(repair, bool):
        raise TypeError("repair must be a boolean")
    std: Word = []
    n = 2 * mg

    while st != std:
        std = st.copy()

        # Original rules 1-2. The bound on k2 changes with len(st); do not use a fixed range.
        for k1 in range(2, n + 1):
            my = list(range(k1 - 1, 0, -1)) + list(range(-n, -k1))
            k2 = 1
            while 2 + k2 * (n - 1) <= len(st):
                fr = [k1] + my * k2 + [-k1]
                # The relation A*a_k*B = Reverse[B]*a_k*Reverse[A] gives
                # a_k*(Reverse[A]*B^(-1))*a_k^(-1)
                # = Reverse[B]^(-1)*A = Reverse[my].
                # Thus the k2-th power of my must map to the k2-th power of Reverse[my].
                to = (my[::-1] if repair else my) * k2
                st = ListPartReplace(st, fr, to)
                st = ListPartReplace(st, _inverse(fr), _inverse(to))
                k2 += 1

        # Original rules 3-4.
        fr = list(range(-n, 0))
        to = fr[::-1]
        st = ListPartReplace(st, fr, to)
        st = ListPartReplace(st, [-x for x in to], [-x for x in fr])

        # Original rules 5-6.
        for k1 in range(2, n + 1):
            fr = list(range(-k1, -n - 1, -1)) + list(range(1, k1))
            to = fr[::-1]
            st = ListPartReplace(st, fr, to)
            st = ListPartReplace(st, [-x for x in to], [-x for x in fr])

        # Free cancellation, preserving the original order by generator index.
        for k1 in range(1, n + 1):
            st = ListPartReplace(st, [k1, -k1], [])
            st = ListPartReplace(st, [-k1, k1], [])

    return st


def GetConjugator(
    wd: Iterable[int], mg: int, *, repair: bool = True
) -> List[Word]:
    """Return [mn, mz]: a representative and a conjugator, with mn = mz * wd * mz^(-1).

    Follow the original method of extracting a core word from the reduced
    words for wd, wd^2, and wd^3. The conjugator mz need not be shortest.

    repair=True corrects the reduction target and the omitted exponent in
    the special branch. repair=False preserves the original reduction rules
    and power branch solely for comparison with the original code; because
    those rules are flawed, this mode need not satisfy the conjugacy equation.
    Both modes fix the extra list nesting in the original final line.
    See the module documentation for details.
    """
    wd, mg = _checked_word(wd, mg)
    if not isinstance(repair, bool):
        raise TypeError("repair must be a boolean")
    st = MyReduction(wd, mg, repair=repair)
    if not st:
        return [[], []]

    st2 = MyReduction(st * 2, mg, repair=repair)
    st3 = MyReduction(st * 3, mg, repair=repair)

    # st = xl + xr; remove prefix and suffix lengths shared by st2 and st3.
    xl_length = _common_prefix_length(st, st2)
    xr = st[xl_length:]
    # Use explicit end indices: Python [:-0] would incorrectly give an empty list.
    xp = st2[xl_length:len(st2) - len(xr)]
    xpp = st3[xl_length:len(st3) - len(xr)]

    xpl_length = _common_prefix_length(xp, xpp)
    xpr = xp[xpl_length:]
    mw = xpp[xpl_length:len(xpp) - len(xp) + xpl_length]
    mz = xpr + xr

    # Select the last cyclic rotation in the original ordering and track the conjugator.
    mn, shift = _max_rotation(mw)
    mz = _inverse(mw[:shift]) + mz

    block_length = 2 * mg - 1
    if len(mn) % block_length:
        return [mn, mz]
    tt = mn[:block_length]
    power = len(mn) // block_length
    if tt * power != mn:
        return [mn, mz]

    # Preserve the original t and tb construction: tb uses only the first 4*mg rows of t.
    generators = list(range(1, 2 * mg + 1))
    cycle = [-x for x in generators] + generators
    t = [_rotate_left(cycle, i) for i in range(1, 4 * mg + 1)]
    t += [_inverse(row) for row in t]
    tb = [row[:block_length] for row in t[:4 * mg]]

    match: Optional[Tuple[int, int, Word]] = None
    for k1 in range(1, block_length + 1):
        rotated = _rotate_left(tt, k1)
        if rotated in tb:
            match = (k1, tb.index(rotated), rotated)
            break
    if match is None:
        return [mn, mz]

    k1, row_index, rotated = match
    reversed_block = rotated[::-1]
    mn2, shift2 = _max_rotation(reversed_block)
    candidate = mn2 * power if repair else mn2

    # Mathematica Sort compares lists by length first, then by their entries.
    # The candidates have equal length by default; repair=False preserves the length ordering.
    if (len(mn), mn) >= (len(candidate), candidate):
        return [mn, mz]

    # The original Join[{list, integer}, ...] incorrectly introduced a nested list.
    # Concatenate the inverse word, a single-letter list, and the remaining conjugator.
    conjugator = (
        _inverse(reversed_block[:shift2])
        + [-t[row_index][2 * mg - 1]]
        + _inverse(tt[:k1])
        + mz
    )
    return [candidate, MyReduction(conjugator, mg, repair=repair)]


# Also provide aliases following Python naming conventions.
list_power = ListPower
list_part_replace = ListPartReplace
my_reduction = MyReduction
get_conjugator = GetConjugator


def main() -> None:
    # Edit wd and mg here to run your own example.
    wd = [1, 2, 3, 4, 6, -1, -2, -4]
    mg = 3
    print("GetConjugator:", GetConjugator(wd, mg))
    print("MyReduction:", MyReduction(wd, mg))


if __name__ == "__main__":
    main()
