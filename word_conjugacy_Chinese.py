#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""word-conjugacy.nb 的 Python 版本（Python 3.8+，仅使用标准库）。

字用整数列表表示：i 表示第 i 个生成元，-i 表示其逆，[] 表示单位元。
mg 沿用原程序，字母范围为 ±1, ..., ±(2*mg)。
约化规则对应的群关系为 a1*a2*...*a(2g) = a(2g)*...*a2*a1。

保留四个原函数名：ListPower、ListPartReplace、MyReduction、GetConjugator。
输入列表不会被修改。可以直接运行本文件，也可以从其他 Python 文件导入。

用法：
    from word_conjugacy import GetConjugator, MyReduction
    wd = [1, 2, 3, 4, 6, -1, -2, -4]
    mn, mz = GetConjugator(wd, 3)
    # mn：算法选出的共轭代表字；mz：共轭子，约定 mn = mz * wd * mz^(-1)。

转换说明：
1. 保留原程序的替换顺序；一次替换后扫描位置只向右移动一位。
2. 原文 Sort[...][[-1]] 取排序最后一项，因此轮换取数值字典序最大者；
   多个轮换相同时保留最大的轮换步数。这里没有按注释改成 min。
3. 默认 repair=True，修正以下三处问题：
   (a) MyReduction 的规则 1–2：to 改为 Reverse[my]^k2，而不是 my^k2。
       这与其他规则使用的关系 a1...a(2g) = a(2g)...a1 一致。
   (b) GetConjugator 最后一个 Join 中多余的列表嵌套改为一维连接。
   (c) 特殊分支处理 tt^q 时，反向候选字也保留 q 次幂。
4. MyReduction(..., repair=False) 保留原约化规则；
   GetConjugator(..., repair=False) 保留原约化规则及原幂次分支行为，
   但仍修正 (b)，以保证返回一维整数列表。这一模式用于对照原笔记本，
   不应依赖其输出做群元素相等或共轭判断。
   原程序比较不同长度的列表时先比较长度，再比较各项；兼容模式保留此规则。
5. 原规则问题的具体例子（mg=2）：
   MyReduction([1,2,3,4,-1,-2,-3,-4], 2, repair=False)
   返回 [1,3,4,-1,-3,-4]；这个输入是上述关系的关系字，应当代表单位元。
   默认修正版对此输入返回 []。
6. MyReduction 和 GetConjugator 仍使用原笔记本的整体算法；程序转换及测试
   不是对该算法在所有输入上给出唯一共轭正规形式的数学证明。

排序语义参考：
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
        raise TypeError("{} 必须是整数".format(name))
    return int(value)


def _integers(values: Iterable[int], name: str) -> Word:
    return [_integer(x, name + " 中的元素") for x in values]


def _checked_word(values: Iterable[int], mg: int) -> Tuple[Word, int]:
    mg = _integer(mg, "mg")
    if mg < 1:
        raise ValueError("mg 必须是正整数")
    word = _integers(values, "字")
    if any(x == 0 or abs(x) > 2 * mg for x in word):
        raise ValueError("字母必须属于 ±1, ..., ±{}".format(2 * mg))
    return word, mg


def _inverse(word: Sequence[int]) -> Word:
    """字的逆：逆序并逐项取相反数。"""
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
    """对应原 Sort[Table[Append[RotateLeft[...], {i}], ...]][[-1]]。"""
    if not word:
        raise ValueError("原算法没有提取到非空核心字，无法选择轮换")
    # 步数为 1,...,len(word)，不是 0,...,len(word)-1。
    # tuple 的第二项用于复现原程序相同轮换的排序规则。
    return max((_rotate_left(word, i), i) for i in range(1, len(word) + 1))


def ListPower(lt: Iterable[int], pr: int) -> Word:
    """对应 ListPower[lt, pr]：连接重复；负指数取逆；零次幂为 []。"""
    lt = _integers(lt, "lt")
    pr = _integer(pr, "pr")
    nl = lt * abs(pr)
    return _inverse(nl) if pr < 0 else nl


def ListPartReplace(li: Iterable[int], fr: Iterable[int], to: Iterable[int]) -> Word:
    """对应 ListPartReplace[li, fr, to]，从左到右扫描并替换连续子字。

    与原 For 循环一致：替换完成后索引增加 1，并按替换后的长度判断终止。
    本函数只做一轮扫描；MyReduction 的外层循环负责重复约化。
    fr 不能为空。用于约化的替换都不会增加字长。
    """
    st, fr, to = list(li), list(fr), list(to)
    if not fr:
        raise ValueError("fr 不能为空列表")
    k1 = 0
    while k1 <= len(st) - len(fr):
        if st[k1:k1 + len(fr)] == fr:
            st[k1:k1 + len(fr)] = to
        k1 += 1
    return st


def MyReduction(li: Iterable[int], mg: int, *, repair: bool = True) -> Word:
    """对应 MyReduction[li, mg]，按原笔记本的规则顺序反复约化。

    repair=True 修正规则 1–2 的目标字；False 保留原规则，供逐项对照。
    """
    st, mg = _checked_word(li, mg)
    if not isinstance(repair, bool):
        raise TypeError("repair 必须是布尔值")
    std: Word = []
    n = 2 * mg

    while st != std:
        std = st.copy()

        # 原规则 1–2。k2 的上界随着 st 的长度改变，不能换成固定 range。
        for k1 in range(2, n + 1):
            my = list(range(k1 - 1, 0, -1)) + list(range(-n, -k1))
            k2 = 1
            while 2 + k2 * (n - 1) <= len(st):
                fr = [k1] + my * k2 + [-k1]
                # 由 A*a_k*B = Reverse[B]*a_k*Reverse[A] 得
                # a_k*(Reverse[A]*B^(-1))*a_k^(-1)
                # = Reverse[B]^(-1)*A = Reverse[my]。
                # 因而 my 的 k2 次幂也必须替换为 Reverse[my] 的 k2 次幂。
                to = (my[::-1] if repair else my) * k2
                st = ListPartReplace(st, fr, to)
                st = ListPartReplace(st, _inverse(fr), _inverse(to))
                k2 += 1

        # 原规则 3–4。
        fr = list(range(-n, 0))
        to = fr[::-1]
        st = ListPartReplace(st, fr, to)
        st = ListPartReplace(st, [-x for x in to], [-x for x in fr])

        # 原规则 5–6。
        for k1 in range(2, n + 1):
            fr = list(range(-k1, -n - 1, -1)) + list(range(1, k1))
            to = fr[::-1]
            st = ListPartReplace(st, fr, to)
            st = ListPartReplace(st, [-x for x in to], [-x for x in fr])

        # 自由消去，保持原程序按生成元编号逐个处理的顺序。
        for k1 in range(1, n + 1):
            st = ListPartReplace(st, [k1, -k1], [])
            st = ListPartReplace(st, [-k1, k1], [])

    return st


def GetConjugator(
    wd: Iterable[int], mg: int, *, repair: bool = True
) -> List[Word]:
    """返回 [mn, mz]，分别为代表字和共轭子，mn = mz * wd * mz^(-1)。

    使用原程序从 wd、wd^2、wd^3 的约化字中提取核心字的方法。
    这里不保证 mz 是最短共轭子。

    repair=True 修正约化目标字和特殊分支的幂次遗漏。
    repair=False 保留原约化规则和幂次分支，仅供对照原代码；
    由于原规则存在问题，这一模式的结果不保证满足上述共轭等式。
    两种模式都修正原程序最后一行的列表嵌套，详见文件开头的说明。
    """
    wd, mg = _checked_word(wd, mg)
    if not isinstance(repair, bool):
        raise TypeError("repair 必须是布尔值")
    st = MyReduction(wd, mg, repair=repair)
    if not st:
        return [[], []]

    st2 = MyReduction(st * 2, mg, repair=repair)
    st3 = MyReduction(st * 3, mg, repair=repair)

    # st = xl + xr；st2 和 st3 去掉相同长度的首尾部分。
    xl_length = _common_prefix_length(st, st2)
    xr = st[xl_length:]
    # 使用显式右端点，避免 Python 的 [:-0] 意外返回空列表。
    xp = st2[xl_length:len(st2) - len(xr)]
    xpp = st3[xl_length:len(st3) - len(xr)]

    xpl_length = _common_prefix_length(xp, xpp)
    xpr = xp[xpl_length:]
    mw = xpp[xpl_length:len(xpp) - len(xp) + xpl_length]
    mz = xpr + xr

    # 所有循环轮换中取原程序排序的最后一项，同时记录共轭子。
    mn, shift = _max_rotation(mw)
    mz = _inverse(mw[:shift]) + mz

    block_length = 2 * mg - 1
    if len(mn) % block_length:
        return [mn, mz]
    tt = mn[:block_length]
    power = len(mn) // block_length
    if tt * power != mn:
        return [mn, mz]

    # 保留原 t、tb 的构造。原 tb 只取 t 的前 4*mg 行。
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

    # Mathematica 的 Sort 比较列表时先比较长度，再比较元素。
    # 默认两候选等长；repair=False 时也保持原来的不同长度排序。
    if (len(mn), mn) >= (len(candidate), candidate):
        return [mn, mz]

    # 原 Join[{列表, 整数}, ...] 错误地产生了嵌套列表。
    # 这里直接连接逆字、单字母列表及余下的共轭子。
    conjugator = (
        _inverse(reversed_block[:shift2])
        + [-t[row_index][2 * mg - 1]]
        + _inverse(tt[:k1])
        + mz
    )
    return [candidate, MyReduction(conjugator, mg, repair=repair)]


# 同时提供符合 Python 命名习惯的别名。
list_power = ListPower
list_part_replace = ListPartReplace
my_reduction = MyReduction
get_conjugator = GetConjugator


def main() -> None:
    # 修改这里的 wd、mg 即可计算自己的例子。
    wd = [1, 2, 3, 4, 6, -1, -2, -4]
    mg = 3
    print("GetConjugator:", GetConjugator(wd, mg))
    print("MyReduction:", MyReduction(wd, mg))


if __name__ == "__main__":
    main()
