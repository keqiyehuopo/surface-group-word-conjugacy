# Word Reduction and Conjugacy in Surface Groups

This repository contains supplementary program files for [*Word Length Formulae, Normal Forms, Conjugation and Root-finding Algorithms in Surface Groups*](https://arxiv.org/abs/2511.12862v3), by Ke Wang, Qiang Zhang, and Xuezhi Zhao (arXiv:2511.12862v3).

Mathematica and Python implementations of word reduction and conjugacy computations in closed orientable surface groups, using the symmetric presentation

$$
G_g = \left\langle a_1,\ldots,a_{2g} \middle|
a_1\cdots a_{2g}=a_{2g}\cdots a_1\right\rangle,
\qquad g\ge 2.
$$

The programs compute normal forms of words and select conjugacy-class representatives together with explicit conjugators.

## Files

| File | Role |
| --- | --- |
| [word-conjugacy_Chinese.nb](word-conjugacy_Chinese.nb) | Mathematica notebook retaining the original Chinese comments. |
| [word-conjugacy_English.nb](word-conjugacy_English.nb) | The same Mathematica notebook, with the Chinese comments translated into English. |
| [word_conjugacy_Chinese.py](word_conjugacy_Chinese.py) | Python implementation with Chinese comments and docstrings. |
| [word_conjugacy_English.py](word_conjugacy_English.py) | Python implementation with English comments and docstrings. |

Within each programming language, both versions use the same computational logic. Either Python file can be used independently; the notebooks do not require Python.

## Word representation

The parameter `mg` denotes the genus. A word is a list of integers: `i` represents $a_i$, and `-i` represents $a_i^{-1}$, with $1\le i\le 2g$. The empty list represents the identity: `{}` in Mathematica and `[]` in Python. Inverting a word reverses the list and negates every entry.

`GetConjugator` returns two lists, `mn` and `mz`, satisfying the convention

$$
mn = mz\,w\,mz^{-1}
$$

in the group. Here `mn` is the selected representative and `mz` is a conjugator; `mz` is not required to be a shortest conjugator.

## Length-lexicographical order

The **length-lexicographical order (shortlex)** compares words first by the number of letters: shorter words are smaller. Words of equal length are compared from left to right at the first differing letter. Following [Definition 2.1 of the paper](https://arxiv.org/html/2511.12862v3#S2.SS1), the alphabet order, written from smallest to largest in this README's notation, is

$$
a_{2g}\prec\cdots\prec a_2\prec a_1\prec a_1^{-1}\prec a_2^{-1}\prec\cdots\prec a_{2g}^{-1}.
$$

In the integer encoding, this is `2g ≺ ... ≺ 2 ≺ 1 ≺ -1 ≺ -2 ≺ ... ≺ -2g`. For example, when `g = 2`, `[4] ≺ [4, 4]` by length, and `[4, 1] ≺ [3, 4]` by the first letter.

The normal form is the smallest word representing the given element or conjugacy class in this order. Since the alphabet order reverses ordinary numeric order, selecting the numerically lexicographically largest equal-length list with `Sort[...][[-1]]` selects the smallest word in the paper's order.

## Algorithm

1. `MyReduction` repeatedly applies the rewriting rules from the symmetric presentation, together with free cancellation, until the word stops changing.
2. `GetConjugator` reduces $w$, $w^2$, and $w^3$, then compares their prefixes to extract a core word and track a conjugator.
3. It selects a cyclic rotation using the original notebook's ordering (`Sort[...][[-1]]`) and updates the conjugator accordingly.
4. A special branch handles powers of certain blocks of length $2g-1$, comparing the corresponding reversed-block candidate while preserving the full exponent.

## Functions

| Function | Purpose |
| --- | --- |
| `ListPower(lt, pr)` | Concatenate copies of a word; invert for negative powers and return the empty list for exponent zero. |
| `ListPartReplace(li, fr, to)` | Replace contiguous subwords during one left-to-right scan. |
| `MyReduction(li, mg)` | Compute the reduced normal form by repeated rewriting. |
| `GetConjugator(wd, mg)` | Return the selected representative and a conjugator. |

The table uses Python syntax. Mathematica calls use square brackets.

## Quick start

### Mathematica

Open either `.nb` file and evaluate the function definitions in the first cell. Then run:

```wolfram
w = {1, 2, 3, 4, 6, -1, -2, -4};
MyReduction[w, 3]
GetConjugator[w, 3]
```

The second call returns:

```wolfram
{{6, -1, -2, -4, 1, 2, 3, 4}, {-4, -3, -2, -1}}
```

### Python

Requires Python 3.8 or later and only the standard library. Run the built-in example:

```bash
python word_conjugacy_en.py
```

Or import the functions from a file in the same directory:

```python
from word_conjugacy_en import GetConjugator, MyReduction

w = [1, 2, 3, 4, 6, -1, -2, -4]
normal = MyReduction(w, 3)
mn, mz = GetConjugator(w, 3)
print(mn)  # [6, -1, -2, -4, 1, 2, 3, 4]
print(mz)  # [-4, -3, -2, -1]
```

To use the Chinese version, import from `word_conjugacy` instead.

## Implementation notes

Python defaults to `repair=True`. Its optional `repair=False` mode retains the original reduction and power-branch behavior for comparison, while still fixing the nested list. Use the default mode for computations.
