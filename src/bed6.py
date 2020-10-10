#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
后修改的Bed类，满足自己的需求.

2017.3.6 进一步进行修改和完善，补充了正负链的判断过程
2017.4.20 完成了对bed类的修改
2017.6.14 完成对脚本的拆分
共计两个Bed类

Bed6，储存每行Bed文件的类，支持bed3,6,12 以及非正式的4（位置信息，带正负链）
（1）包含正常bed6的各种setters和getters
        set_chrom (get_chrom): 设置（获取）染色体
        set_start (get_start): 设置（获取）起始位点
        set_end (get_end): 设置（获取）终止位点
        set_name (get_name): 设置（获取）名称信息
        set_score (get_score): 设置（获取）分值信息（该位置可以换成其他任意信息）
        set_strand (get_strand): 设置（获取）链信息
（2）重载了多个神奇方法
        __gt__、__eq__、__lt__、__ge__、__le__ 用来根据位置信息（染色体，起始终止位点，和正负链）判断Bed6的大小
        __hash__ 用来获取hash值，因此Bed6是iterable，支持in和not in（没有单独在重载__contain__等）等
        __and__ 用来判断Bed6之间是否具有重合位点
        __mul__ 用来判断Bed6的范围是否完全覆盖第二个Bed6
        __mod__ 用来计算两个Bed6的start site的距离
        __truediv__ 用来计算两个Bed6的中心位点的距离
（3）添加新方法
        cal_dist: 参数center（默认为False）, 用来根据参数，自动计算Bed6之间相应的距离
        __generate_chr_order: 自动生成chr1~chr22的染色体list，其他多余的chr会按照字符串的顺序导入
                              主要用于支持大小对比，并且用于Bed6的sort
        __get_index: 通过内置的__chrom来获取在染色体list的位置
"""
import re
import sys
import typing
from tqdm import tqdm

class Bed6():
    u"""Bed class."""

    __chrom = '.'
    __start = -1
    __end = -1
    __name = '.'
    __score = '.'
    __strand = '.'
    __append = ''  # 除了标准六行之外的信息
    __addition_chrom = []
    __feature = ""
    original = ""

    def __init__(self, chrom: str=None, start: int=None, end: int=None, name: str=None, score: str=None,
                 strand: str=None, feature:str=None, append=None, original: str = None):
        u"""Init."""
        self.__order = self.__generate_chr_order()  # 必须先生成__order，随后需要用
        self.set_chrom(chrom)
        self.set_start(start)
        self.set_end(end)
        self.set_name(name)
        self.set_score(score)
        self.set_strand(strand)
        self.__feature = feature
        self.original = original
        if append:
            if isinstance(append, list):
                append = '\t'.join(append)
            self.__append = append

    @staticmethod
    def __generate_chr_order():                     # 生成染色体的顺序
        order = []
        for i in range(1, 23):
            order.append('chr' + str(i))
        order.append('chr' + 'X')
        order.append('chr' + 'Y')
        order.append('chr' + 'M')
        return order

    # 各路setters
    def set_chrom(self, chrom):
        u"""Set chromosome."""
        if chrom is None:
            return
        if chrom not in self.__order:
            if chrom not in self.__addition_chrom:
                # 如果染色体不在默认名单里，就添加到额外的名单里
                self.__addition_chrom.append(chrom)
        self.__chrom = chrom

    def set_start(self, start):
        u"""Set start site."""
        if start is None:
            return
        if isinstance(start, int):
            self.__start = start
        else:
            self.__start = int(start)

    def set_end(self, end):
        u"""Set end site."""
        if end is None:
            return
        # 修改，完成，原始数据先转换成合适的类型，然后进行比对
        if isinstance(end, int):
            self.__end = end
        else:
            self.__end = int(end)
        if self.__end < self.__start:
            raise ValueError('end < start')

    def set_name(self, name):
        u"""Set track name."""
        if name is not None:
            self.__name = name

    def set_score(self, score):
        u"""Set track score, but also can set as other message."""
        if score is not None:   # 应该就是这出的问题
            self.__score = score

    def set_strand(self, strand):
        u"""Set track strand + or - or . ."""
        if strand is not None and strand in ('+', '-', '.'):
            self.__strand = strand

    def set_append(self, append):
        u"""Set append msg."""
        self.__append = append

    # 各种getters
    @property
    def chrom(self):
        u"""Return chromosome."""
        return self.__chrom

    @property
    def start(self):
        u"""Return track start site."""
        return self.__start

    @property
    def end(self):
        u"""Return track end site."""
        return self.__end

    @property
    def name(self):
        u"""Return track name."""
        return self.__name

    @property
    def score(self):
        u"""Return the message at column 5."""
        return self.__score

    @property
    def strand(self):
        u"""Return track strand."""
        return self.__strand

    @property
    def feature(self):
        return self.__feature

    @property
    def append(self):
        u"""Return append msg."""
        return self.__append

    @property
    def is_empty(self) -> bool:
        return str(self) == str(Bed6())

    def get_bed(self):
        u"""Return track."""
        standard = ""
        for i in [self.__chrom, self.__start, self.__end, self.__name,
                  self.__score, self.__strand]:
            if i != '':                                     # 如果有一项为空，就不输出了
                if standard:
                    standard = standard + str(i) + '\t'
                else:
                    standard = str(i) + '\t'
        # 除去多余的\t，经测试，使用了正则才成功去除，rstrip和strip都失败，原因不明，可能我这不算字符串吧
        standard = re.sub(r'\t$', '', standard)
        # 原本standard和append之间还存在一个多余的\t，现在去除掉
        if self.__append:
            standard = standard + '\t' + self.__append
        return standard

    def get_bed4(self) -> str:
        return f"{self.__chrom}\t{self.__start}\t{self.__end}\t{self.__strand}"

    def __str__(self) -> str:
        return self.get_bed()

    def get_center(self):
        u"""返回这个区域的中心."""
        return (self.__end - self.__start) / 2 + self.__start

    # 下面重载几个常用的选项，来满足自己的比较需求
    def __get_index(self, chrom):
        u"""这个get_index与Bed的不一样，是用来比较chrom谁比较靠前的."""
        try:
            if chrom in self.__order:              # 染色体正常就通过这个获取制定的排列顺序
                return self.__order.index(chrom)
            else:
                self.__addition_chrom.sort()       # 不常见染色体，就按字符串排序，意思意思
                self.__order.extend(self.__addition_chrom)
                return self.__order.index(chrom)
        except ValueError:
            print(chrom)
            sys.exit(1)

    def __lt__(self, other):                    # 大于, isupperstream
        u"""定义大于，用于比较两个bed track的先后循序，按照染色体，起始和终止位点， 上游."""
        if self.__get_index(self.__chrom) < self.__get_index(other.__chrom):
            return True
        elif self.__chrom == other.__chrom:
            if self.__end < other.__start:
                return True
        return False

    def __eq__(self, other):                # 相等
        u"""定义两条track相等的情况，这里只判断了链的位置信息，并不判断名称等等."""

        # 2017.8.18记，添加链的判断，和修复end位点的判断
        if self.__chrom == other.__chrom and \
                self.__start == other.__start and \
                self.__end == other.__end and \
                self.__compare_strand__(other):
            return True
        return False

    def __gt__(self, other):                # 小于  isdownstream
        u"""定义小于的情况, 下游."""
        if self.__get_index(self.__chrom) > self.__get_index(other.__chrom):
            return True
        elif self.__chrom == other.__chrom:
            if self.__start > other.__end:
                return True
        return False

    def __le__(self, other):
        u"""
        定义大于等于，用于比较两个bed track的先后循序，按照染色体，起始和终止位点.

        大于等于，两个区域可有可无交集，但是self的起始位点>=other，并且添加链的判断
        """
        if self.__get_index(self.__chrom) < self.__get_index(other.__chrom):
            return True
        elif self.__chrom == other.__chrom:
            if self.__start <= other.__start:
                if self.__compare_strand__(other) or \
                        (self.__strand == '+' or other.__strand == '-'):
                    return True
        return False

    def __ge__(self, other):            # 小于等于，两个区域可有可无交集，但是self的终止位点<=other
        u"""定义小于等于的情况."""
        if self.__get_index(self.__chrom) > self.__get_index(other.__chrom):
            return True
        elif self.__chrom == other.__chrom:
            if self.__end >= other.__end:
                if self.__compare_strand__(other) or \
                        (self.__strand == '-' or other.__strand == '+'):
                    return True
        return False

    def __compare_strand__(self, other):
        u"""用来判断链是否正确."""
        if self.__strand == other.__strand or\
                self.__strand in ('', '.') or\
                other.__strand in ('', '.'):
            return True
        return False

    def __add__(self, other):
        u""" 用+号来merge两个不同的区域 """

        self.__start =  min(self.__start, other.__start)
        self.__end = max(self.__end, other.__end)
        return self

    def __and__(self, other):
        u"""用 & 来判断两个bed之间是否重合."""
        # 必须判断两个区域的染色体和正负链问题
        if self.__chrom == other.__chrom:
            if self.__compare_strand__(other):
                if self.__start <= other.__end and\
                        self.__end >= other.__start:
                    return True
                elif other.__start <= self.__end and\
                        other.__end >= self.__start:
                    return True
        return False

    def __hash__(self):
        u"""只有指定了hash方式，才能用于等于比较，用于in, not in等."""
        return hash(self.get_bed())

    def __mul__(self, other):
        u"""用来统计该bed是否包含有另外一个bed类的范围的."""
        if self.__chrom == other.__chrom and \
            self.__compare_strand__(other) and \
            self.__start <= other.__start and\
                self.__end >= other.__end:
            return True

        return False

    def __mod__(self, other):
        u"""取余 %，用来返回两个start site的距离."""
        return self.__start - other.__start

    def __truediv__(self, other):
        u"""除法/，用来返回两个bed中心的距离."""
        return self.get_center() - other.get_center()

    def isdownstream_same(self, other):
        u"""供downstream_target使用的判断上下游用的."""
        if self.__chrom == other.__chrom:
            if self.__start <= other.__end:
                return True
            else:
                return 2
        elif self.__get_index(self.__chrom) > self.__get_index(other.__chrom):
            return 2
        else:
            return 3

    def is_same_chrom(self, other):
        u"""
        用来判断两个bed6的染色体的大小
        """
        if self.__chrom == other.__chrom:
            return True
        elif self.__get_index(self.__chrom) > self.__get_index(other.__chrom):
            return 1
        else:
            return -1

    def overlap(self, other):
        dist = -1
        if self.__chrom == other.__chrom:
            dist = min(self.__end, other.__end) - max(self.__start, other.__start)
            if dist <= 0:
                dist = -1
        return dist


class BedPair(object):
    def __init__(self, src: Bed6, target: Bed6=None):
        self.src = src
        self.target = target if target is  not None else Bed6()

    def __str__(self) -> str:
        dist = self.src.overlap(self.target)
        if dist < 0:
            dist = 0
        
        target = str(self.target) if self.target.is_empty else self.target.original
        return f"{self.src.original}\t{target}\t{dist}"

    def __gt__(self, other):
        if self.src != other.src:
            return self.src > other.src
        return self.target > other.target

    def __lt__(self, other):
        if self.src != other.src:
            return self.src < other.src
        return self.target < other.target

    @classmethod
    def save(cls, data, output: str, keep_all=False):
        u""" save to file """
        with open(output, "w+") as w:
            for i in sorted(data):
                if not keep_all and i.target.is_empty:
                    continue
                w.write(str(i) + "\n")


class BedUtils(object):
    def __init__(self) -> None:
        u""" Nothing to do """
        pass

    @classmethod
    def read(cls, input_file: str, required_regions: typing.List[Bed6] = None) -> typing.List[Bed6]:
        res = []

        required = {}
        if required_regions:
            for i in required_regions:
                temp = required.get(i.chrom, [0, 0])
                temp[0] = min(temp[0], i.start) if temp[0] else i.start
                temp[1] = max(temp[1], i.end) if temp[1] else i.end
                required[i.chrom] = temp

        with open(input_file) as r:
            for line in tqdm(r, desc="reading " + input_file):
                if line.startswith("#") or line.startswith("track"):
                    continue
                
                lines = line.strip().split("\t")

                if required:
                    if not lines[0] in required.keys():
                        continue
                        
                    if int(lines[2]) < required[lines[0]][0]:
                        continue
                        
                    if int(lines[1]) > required[lines[0]][1]:
                        continue

                name, score, strand, append = None, None, None, None
           
                if len(lines) == 4:
                    if lines[3] in ["+", "-"]:
                        strand = lines[3]
                    else:
                        score = lines[3]
                elif len(lines) == 5:
                    name, score = lines[3], lines[4]
                elif len(lines) >= 6:
                    name, score, strand = lines[3], lines[4], lines[5]
                    append = "\t".join(lines[6:])

                res.append(
                    Bed6(
                        chrom=str(lines[0]), 
                        start=int(lines[1]),
                        end=int(lines[2]),
                        name=name,
                        score=score,
                        strand=strand,
                        append=append,
                        original=line.strip()
                    )
                )

        return res
    
    @classmethod
    def write(cls, data, output: str):
        with open(output, "w+") as w:
            for i in sorted(data):
                w.write(str(i) + "\n")

    @classmethod
    def __list2dict__(cls, data:typing.List[Bed6], res: typing.Dict[str, typing.List[Bed6]]=None, strandness: bool = False) -> typing.Dict[str, typing.List[Bed6]]:
        if res is None:
            res = {}
        
        for i in data:
            key = i.chrom if strandness else "%s#%s" % (i.chrom, i.strand)
            temp = res.get(key, [])
            temp.append(i)
            res[key] = temp
        
        return res

    @classmethod
    def is_upstream(cls, src: Bed6, target: Bed6, distance: int=0) -> bool:
        u"""
        determine whether src is upstream of target
        """
        if src.chrom != target.chrom:
            return src.chrom < target.chrom
        
        return src.end + distance < target.start
    
    @classmethod
    def is_downstream(cls, src: Bed6, target: Bed6, distance: int=0) -> bool:
        u""" determine whether src is downstream of target """
        if src.chrom != target.chrom:
            return src.chrom > target.chrom
        
        return src.start > target.end + distance

    @classmethod
    def is_cover(cls, src: Bed6, target: Bed6, strandness: bool=True) -> bool:
        res = src.chrom != target.chrom and src.start <= target.start and src.end >= target.end

        if not strandness:
            res = res and src.strand == target.end

        return res

    @classmethod
    def sort(cls, beds: typing.List[Bed6]) -> typing.List[Bed6]:
        u""" as name says """
        return sorted(beds, key=lambda x:[x.chrom, x.start, x.end, x.strand])

    @classmethod
    def self_merge(cls, beds: typing.List[Bed6], distance: int=0, strandness: bool = True) -> typing.List[Bed6]:
        u"""
        merge a list of beds 
        """
        if len(beds) < 2:
            return beds

        data = cls.__list2dict__(beds, strandness=strandness)

        res = []
        for vals in data.values():
            vals = cls.sort(vals)
            last = vals[0]

            for i in range(1, len(vals)):
                if cls.is_upstream(last, vals[i], distance=distance):
                    res.append(last)
                    last = vals[i]
                elif cls.is_downstream(last, vals[i], distance=distance): #  distance=distance
                    res.append(vals[i])
                else:
                    last += vals[i]
            
            res.append(last)
        
        return cls.sort(res)

    @classmethod
    def __intersect__(cls, src: Bed6, target: Bed6) -> typing.Optional[Bed6]:
        u""" extract intersects from two regions """

        if not src & target:
            return None

        return Bed6(
            chrom=src.chrom, feature=src.feature,
            start=max(src.start, target.start),
            end=min(src.end,  target.end),
            strand=src.strand
        )

    @classmethod
    def intersect(cls, src: typing.List[Bed6], target: typing.List[Bed6], strandness: bool=True) -> typing.List[Bed6]:
        u"""  """
        res = []

        ref = cls.__list2dict__(src, strandness=strandness)
        data = cls.__list2dict__(target, strandness=strandness)
        
        for key in set(ref.keys()) & set(data.keys()):
            ref_lst = cls.sort(ref[key])
            data_lst = cls.sort(data[key])

            i, j = 0, 0
            first_match = -1
            while i < len(ref_lst) and j < len(data_lst):
                if cls.is_upstream(ref_lst[i], data_lst[j]):
                    i += 1

                    if first_match >= 0:
                        j = first_match
                        first_match = -1
                elif cls.is_downstream(ref_lst[i], data_lst[j]):
                    j += 1
                else:
                    if first_match < 0:
                        first_match = j

                    temp = cls.__intersect__(ref_lst[i], data_lst[j])
                    if temp is not None:
                        res.append(temp)
                    j +=1

        return cls.sort(res)

    @classmethod
    def __substract__(cls, src: Bed6, target: Bed6) -> Bed6:
        u""" substract unique part """
        return Bed6(
            chrom=src.chrom, feature=src.feature,
            start=max(src.start, target.start),
            end=min(src.end, target.end),
            strand=src.strand
        )

    @classmethod
    def subtract(cls, src: typing.List[Bed6], target: typing.List[Bed6], strandness: bool=True) -> typing.List[Bed6]:
        u"""  """
        res = []
        ref = cls.__list2dict__(src, strandness=strandness)
        data = cls.__list2dict__(target, strandness=strandness)

        for key in set(ref.keys()) & set(data.keys()):
            ref_lst = cls.self_merge(ref[key])
            data_lst = cls.self_merge(data[key])

            i, j = 0, 0
            first_match = -1
            overlapped = set()
            while i < len(ref_lst) and j < len(data_lst):
                if cls.is_upstream(ref_lst[i], data_lst[j]):

                    if ref_lst[i] not in overlapped:
                        res.append(ref_lst[i])
         
                    if first_match >= 0:
                        j = first_match
                        first_match = -1
                    i += 1
                elif cls.is_downstream(ref_lst[i], data_lst[j]):
                    j += 1
                else:
                    if first_match < 0:
                        first_match = j
                    
                    overlapped.add(ref_lst[i])
                    
                    if cls.is_cover(ref_lst[i], data_lst[j], strandness=strandness):
                        res.append(cls.__substract__(ref_lst[i], data_lst[j]))
                    j +=1

        return res

    @classmethod
    def overlap(cls, src: typing.List[Bed6], target: typing.List[Bed6], strandness: bool=True) -> typing.List[BedPair]:
        u"""
        find the overlapped region and output it
        """
        res = []
        src_collection, target_collection = {}, {}
        for i in src:
            key = f"{i.chrom}{i.strand}" if strandness else i.chrom
            temp = src_collection.get(key, [])
            temp.append(i)
            src_collection[key] = temp

        for i in target:
            key = f"{i.chrom}{i.strand}" if strandness else i.chrom
            temp = target_collection.get(key, [])
            temp.append(i)
            target_collection[key] = temp

        matched = set()
        for key in set(src_collection.keys()) & set(target_collection.keys()):
            temp_src = sorted(src_collection[key])
            temp_target = sorted(target_collection[key])

            i, j, first = 0, 0, None
            while  i < len(temp_src) and j < len(temp_target):
                if temp_src[i] > temp_target[j]:
                    j += 1
                elif temp_src[i] < temp_target[j]:
                    if temp_src[i] not in matched:
                        res.append(BedPair(temp_src[i]))
                    if first is not None:
                        j = first
                    i += 1
                else:
                    if first is None:
                        first = j
                    if temp_src[i].overlap(temp_target[j]) >= 0:
                        res.append(BedPair(temp_src[i], temp_target[j]))
                        matched.add(temp_src[i])

                    j += 1
        return res


if __name__ == '__main__':
    i = "/mnt/raid61/Personal_data/zhangyiming/code/afe/bulk/mountainClimber/tests/climb/S0168J11_EB_RU.bed"
    t = "/mnt/raid61/Personal_data/zhangyiming/code/afe/bulk/mountainClimber/tests/climb/tus_merged.annot.Homo_sapiens.GRCh38.93_singleGenes.bed"

    i = BedUtils.read(i)
    t = BedUtils.read(t)
    print(len(i))
    print(len(t))
    BedPair.save(BedUtils.overlap(i, t, strandness=False), "/mnt/raid61/Personal_data/zhangyiming/code/afe/bulk/mountainClimber/tests/bed_wo.sample")


