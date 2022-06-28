import logging
import os
from collections import defaultdict
from collections import namedtuple

import pandas as pd

RESULTSPATH = "./results"
COLBARCODENAME = "Barcode"
COLINFONAME = "Info"
# debug
DEBUG = True
DEBUGFILE = "test.csv"


# TODO: Recursive file search functionality
# that only picks final_indel_results.tsv

class LengthAnalyzer():
    def __init__(self):
        if DEBUG:
            self.path = os.path.join(RESULTSPATH, DEBUGFILE)
        else:
            self.path = RESULTSPATH

        self.storage = defaultdict(lambda: "Barcode Error")
        self.ParsedIndelInfo = namedtuple("ParsedIndelInfo",
                                          ["info", "count", "ratio"])  # ratio = certain indel count / total indel count
        self.FlabbyIndelInfo = namedtuple("FlabbyIndelInfo", ["position", "length", "type"])

    def __unpack(self, compact: str):
        """
        :Example: 
        22M1D : mutation on the position 22, and that is the 1 bp insertion

        returns (22, 1, "D")

        :param compact: 
        :type compact: 
        :return: 
        :rtype: 
        """

        # String processing
        position, detail = compact.split("M")
        length, type = detail[:-1], detail[-1]

        position = int(position.strip())
        length = int(length)

        return position, length, type

    def __info_tokenizer(self, info):
        """

        :param info:
        gets a string from a certain barcode 
        :type info: 
        str
        :return:
        Returns a list of named tuple that contains indel patterns 
        :rtype: 
        """
        temp = info.split(',')
        temp = [token.strip() for token in temp]

        parsed = list()
        for item in temp:
            # The Last element
            if len(item) == 0:
                continue

            compact_info, count, ratio = item.split(':')
            flabby_info = self.FlabbyIndelInfo(*self.__unpack(compact_info))
            parsed.append(self.ParsedIndelInfo(flabby_info, float(count), float(ratio)))

        return parsed

    def parse_info(self):
        """
        tsv to dictionary: for in-memory processing
        """

        # TODO: for multiple files

        if DEBUG:
            # Single file version
            df = pd.read_csv(self.path, sep='\t')

        # current progress  
        self.total_length = df.shape[0]

        for index, row in df.iterrows():
            id, info = row[COLBARCODENAME], row[COLINFONAME]
            self.storage[id] = self.__info_tokenizer(info)

            # Displaying current progress
            if index % 1000 == 0:
                print(f"Loading Indel Information... {index} out of {self.total_length}")

    def process(self):

        # 정보가 21M1I:1880:34.0, 17M13D:714:12.0, ... 이런식으로 나열되어 있으니
        # 각 barcode (target sequence)의 info 내에서
        # (1) M~I:xxxx 와 M~D:yyyy 인 것으로 나누고 x와 y 값의 합 구하기
        # (2) M~I:xxxx 들을 모두 모아놓고... ~ (insertion 길이) x (해당 길이 insertion을 가진 read의 비율; xxxx/(xxxx summation))
        # (3) (2)와 동일하지만, deletion으로

        logging.basicConfig(format="%(messages)s")
        log = logging.getLogger()
        # Loading on the memory -> in memory processing
        # TODO: multiprocessing
        self.parse_info()
        print("Data Load Completed")

        # Instantiate useful variables
        # (1) insertion, deletion의 비율 dictionary
        ins_ratios = dict.fromkeys(self.storage.keys())
        del_ratios = dict.fromkeys(self.storage.keys())
        # (2) insertion 평균길이 dictionary; weighted length
        ins_lengths = dict.fromkeys(self.storage.keys())
        # (3) deletion의 평균길이 dictionary; weighted length
        del_lengths = dict.fromkeys(self.storage.keys())

        mutation_collection = list()
        for cur, barcode in enumerate(self.storage):
            # The core processing part

            insertions = list()
            deletions = list()

            for item in self.storage[barcode]:
                info = getattr(item, "info")
                mt_type = getattr(info, "type")

                assert mt_type == "I" or "D", f"Invalid mutation type: {mt_type}"

                if mt_type == "I":  # Insertions
                    insertions.append(item)
                elif mt_type == "D":  # Deletions
                    deletions.append(item)

            ins_ratios[barcode] = len(insertions) / (len(insertions) + len(deletions))
            del_ratios[barcode] = 1 - ins_ratios[
                barcode]  # Repeated floating point calculations can hamper the performance 

            # Weighted Insertion Length
            weighted_ins_len = 0.0
            weighted_del_len = 0.0
            ratio_guarantee = 0.0

            # TODO: Modularizing
            for occ in insertions:
                weighted_ins_len += getattr(getattr(occ, "info"), "length") * getattr(occ, "ratio")

                # debug

                # ratio_guarantee += getattr(occ, "ratio")

            ins_lengths[barcode] = weighted_ins_len / 100.0

            # Weighted Deletion Length
            for occ in deletions:
                weighted_del_len += getattr(getattr(occ, "info"), "length") * getattr(occ, "ratio")

                # debug
                # ratio_guarantee += getattr(occ, "ratio")

            del_lengths[barcode] = weighted_del_len / 100.0

            # TODO: parametrizing
            # if ratio_guarantee < 100 - 10:
            #     log.warning(f"The sum of indel ratio is not adding up to 100%. Manul review recommended."
            #                 f"Barcode : {barcode}")

            # Displaying current progress
            if cur % 1000 == 0:
                print(f"Loading Indel Information... {cur} out of {len(self.storage)}")

        # TODO: for exporting processed data
        column_names = ["Barcode", "Ins Ratio", "Del Ratio", "Weighted Ins Len", "Weighted Del Len"]

        # Populating a dataframe
        data = [[barcode, ins_ratios[barcode], del_ratios[barcode], ins_lengths[barcode], del_lengths[barcode]] for
                barcode in
                self.storage]
        df = pd.DataFrame(data, columns=column_names)
        print(df)

        df.to_csv(f"{self.path}_analyzed.csv", index=False)


debug = LengthAnalyzer()
# test = debug.info_tokenizer("22M1D:1482:5.6000000000000005, 20M5D:1302:4.9, 16M23D:1168:4.3999999999999995, 24M1D:1121:4.2, 20M10D:1007:3.8, 29M1D:997:3.8, 16M11D:756:2.8000000000000003, 16M12D:688:2.6, 18M15D:649:2.4, 19M11D:648:2.4, 19M6D:634:2.4, 18M25D:541:2.0, 18M12D:498:1.9, 13M23D:352:1.3, 12M15D:347:1.3, 12M32D:347:1.3, 20M9D:346:1.3, 12M22D:319:1.2, 29M17I:308:1.2, 12M12D:292:1.0999999999999999, 18M20D:271:1.0, 8M20D:263:1.0, 18M7D:254:1.0, 25M4D:251:0.8999999999999999, 20M4D:235:0.8999999999999999, 17M14D:227:0.8999999999999999, 11M19D:204:0.8, 15M15D:191:0.7000000000000001, 16M10D:172:0.6, 14M17D:157:0.6, 20M14D:154:0.6, 11M18D:152:0.6, 20M6D:149:0.6, 16M16D:149:0.6, 12M14D:141:0.5, 9M14D:127:0.5, 13M16D:122:0.5, 31M17D:121:0.5, 18M21D:118:0.4, 12M16D:118:0.4, 14M12D:114:0.4, 21M9D:112:0.4, 20M8D:111:0.4, 7M27D:108:0.4, 6M22D:106:0.4, 5M38D:105:0.4, 12M11D:104:0.4, 17M9D:104:0.4, 10M21D:104:0.4, 18M33D:103:0.4, 18M8D:102:0.4, 15M11D:100:0.4, 15M14D:98:0.4, 5M34D:98:0.4, 28M9I:97:0.4, 23M1I:96:0.4, 13M14D:94:0.4, 16M34D:90:0.3, 16M13D:89:0.3, 13M19D:87:0.3, 15M28D:86:0.3, 18M13D:86:0.3, 29M2I:86:0.3, 29M6I:85:0.3, 21M12D:84:0.3, 16M17D:83:0.3, 15M13D:82:0.3, 21M11D:82:0.3, 15M18D:79:0.3, 26M10D:77:0.3, 15M10D:76:0.3, 10M35D:76:0.3, 23M3D:75:0.3, 22M7D:71:0.3, 8M22D:69:0.3, 23M8D:66:0.2, 22M29D:66:0.2, 8M32D:66:0.2, 35M6D:65:0.2, 8M30D:65:0.2, 12M27D:64:0.2, 17M13D:64:0.2, 16M20D:62:0.2, 11M27D:62:0.2, 12M28D:60:0.2, 20M7D:59:0.2, 8M23D:59:0.2, 12M18D:59:0.2, 26M4D:59:0.2, 32M6D:58:0.2, 12M38D:58:0.2, 11M17D:57:0.2, 25M3D:56:0.2, 8M33D:53:0.2, 19M15D:53:0.2, 22M18D:52:0.2, 24M19D:50:0.2, 25M3I:50:0.2, 32M4I:50:0.2, 16M21D:50:0.2, 18M22D:50:0.2, 7M28D:49:0.2, 22M5D:49:0.2, 34M6I:49:0.2, 10M31D:48:0.2, 34M1D:48:0.2, 25M17D:48:0.2, 27M12D:47:0.2, 25M5D:46:0.2, 22M6D:45:0.2, 22M13D:44:0.2, 19M20D:43:0.2, 18M31D:43:0.2, 21M14D:43:0.2, 16M15D:43:0.2, 26M7D:43:0.2, 22M8I:42:0.2, 24M2D:41:0.2, 27M8D:41:0.2, 3M27D:40:0.2, 4M28D:40:0.2, 22M15D:40:0.2, 24M9D:39:0.1, 21M10D:39:0.1, 27M3I:39:0.1, 34M12D:38:0.1, 31M11D:38:0.1, 19M7D:36:0.1, 17M12D:36:0.1, 5M21D:35:0.1, 27M15I:35:0.1, 9M22D:34:0.1, 24M14D:34:0.1, 28M14D:34:0.1, 18M16D:34:0.1, 13M13D:34:0.1, 22M8D:34:0.1, 27M2I:33:0.1, 11M25D:33:0.1, 9M28D:33:0.1, 26M1D:33:0.1, 27M7D:30:0.1, 32M7D:30:0.1, 24M3I:29:0.1, 26M6D:29:0.1, 25M21D:29:0.1, 27M15D:29:0.1, 28M20I:29:0.1, 23M19D:27:0.1, 19M14D:27:0.1, 25M6D:26:0.1, 23M10D:26:0.1, 19M19D:26:0.1, 7M19D:26:0.1, 8M16D:26:0.1, 8M21D:25:0.1, 9M19D:24:0.1, 20M11D:24:0.1, 25M1I:24:0.1, 16M9D:24:0.1, 19M17D:23:0.1, 25M7D:23:0.1, 28M2I:23:0.1, 23M13D:23:0.1, 23M2I:23:0.1, 22M9D:23:0.1, 19M12D:22:0.1, 19M13D:22:0.1, 31M12I:22:0.1, 20M23D:22:0.1, 24M10D:22:0.1, 26M6I:22:0.1, 14M26D:22:0.1, 31M6I:21:0.1, 14M18D:21:0.1, 34M17D:21:0.1, 29M9D:21:0.1, 18M29D:21:0.1, 25M9D:21:0.1, 24M13D:21:0.1, 21M15D:21:0.1, 33M11D:20:0.1, 32M10I:20:0.1, 32M5I:20:0.1, 29M5I:19:0.1, 34M8I:19:0.1, 23M6D:18:0.1, 28M3D:18:0.1, 21M21D:18:0.1, 20M13D:18:0.1, 15M20D:18:0.1, 30M17I:18:0.1, 8M28D:18:0.1, 26M19I:18:0.1, 22M4I:17:0.1, 27M32D:17:0.1, 31M7I:17:0.1, 34M20I:17:0.1, 28M2D:17:0.1, 36M5I:16:0.1, 21M6D:16:0.1, 32M12D:16:0.1, 31M2D:16:0.1, 35M8D:16:0.1, 33M13D:16:0.1, 11M29D:16:0.1, 28M10D:16:0.1, 29M23I:15:0.1, 14M14D:15:0.1, 16M8D:15:0.1, 28M13D:15:0.1, 17M25D:15:0.1, 22M11D:15:0.1, 27M5I:15:0.1, 13M18D:14:0.1, 30M9I:14:0.1, 7M21D:14:0.1, 26M3I:14:0.1, 29M4I:14:0.1, 27M13I:14:0.1, 23M15D:14:0.1, 19M21D:13:0.0, 6M24D:13:0.0, 7M23D:13:0.0, 22M5I:13:0.0, 17M8D:13:0.0, 33M5I:13:0.0, 22M14D:13:0.0, 20M16D:12:0.0, 26M1I:12:0.0, 28M11I:12:0.0, 15M12D:12:0.0, 23M16D:12:0.0, 7M32D:12:0.0, 5M23D:12:0.0, 13M24D:12:0.0, 26M5D:12:0.0, 24M11D:11:0.0, 31M4I:11:0.0, 5M31D:11:0.0, 28M1I:11:0.0, 30M1I:11:0.0, 27M3D:11:0.0, 24M16D:11:0.0, 23M7D:11:0.0, 11M24D:11:0.0, 3M41D:11:0.0, 13M15D:11:0.0, 14M37D:10:0.0, 22M16I:10:0.0, 27M7I:10:0.0, 34M14I:10:0.0, 30M5I:10:0.0, 15M34D:10:0.0, 24M40I:10:0.0, 26M2I:10:0.0, 36M8I:10:0.0, 34M10I:9:0.0, 28M7I:9:0.0, 21M23D:9:0.0, 24M6D:9:0.0, 14M13D:9:0.0, 22M6I:9:0.0, 30M6D:9:0.0, 17M11D:9:0.0, 9M18D:9:0.0, 34M11D:9:0.0, 22M10D:9:0.0, 27M36I:8:0.0, 19M18D:8:0.0, 21M20D:8:0.0, 23M11D:8:0.0, 32M22I:8:0.0, 24M4I:8:0.0, 9M33D:8:0.0, 28M4I:7:0.0, 18M11D:7:0.0, 30M6I:7:0.0, 31M1D:7:0.0, 26M3D:7:0.0, 5M20D:7:0.0, 4M26D:7:0.0, 30M16D:7:0.0, 30M11I:7:0.0, 30M11D:7:0.0, 26M7I:7:0.0, 28M16I:7:0.0, 20M20D:7:0.0, 18M6D:7:0.0, 23M8I:7:0.0, 17M19D:7:0.0, 26M5I:6:0.0, 30M9D:6:0.0, 33M12D:6:0.0, 27M4I:6:0.0, 16M7D:6:0.0, 7M33D:6:0.0, 27M18I:6:0.0, 21M13D:6:0.0, 27M20I:6:0.0, 15M8D:6:0.0, 27M10I:6:0.0, 9M16D:6:0.0, 21M32D:6:0.0, 23M21D:6:0.0, 24M18D:6:0.0, 21M7D:5:0.0, 34M11I:5:0.0, 3M22D:5:0.0, 20M15D:5:0.0, 20M28D:5:0.0, 18M26D:5:0.0, 9M25D:5:0.0, 32M9D:5:0.0, 30M8I:5:0.0, 23M4I:4:0.0, 17M21D:4:0.0, 23M12D:4:0.0, 20M12D:4:0.0, 24M2I:4:0.0, 27M6D:4:0.0, 25M6I:4:0.0, 22M2I:4:0.0, 11M13D:4:0.0, 3M37D:4:0.0, 20M21D:4:0.0, 28M8I:4:0.0, 22M1I:4:0.0, 8M19D:4:0.0, 29M7I:4:0.0, 17M15D:4:0.0, 26M12D:4:0.0, 29M9I:4:0.0, 21M3D:4:0.0, 8M37D:4:0.0, 18M23D:3:0.0, 8M27D:3:0.0, 9M41D:3:0.0, 34M5I:3:0.0, 6M31D:3:0.0, 17M23D:3:0.0, 9M17D:3:0.0, 23M9D:3:0.0, 10M41D:3:0.0, 35M51I:3:0.0, 25M2D:3:0.0, 32M3I:3:0.0, 25M15D:3:0.0, 18M32D:3:0.0, 11M21D:3:0.0, 9M26D:3:0.0, 27M8I:3:0.0, 27M4D:3:0.0, 4M22D:3:0.0, 14M33D:3:0.0, 19M24D:3:0.0, 28M5I:3:0.0, 29M19D:3:0.0, 29M10D:3:0.0, 17M29D:3:0.0, 3M31D:2:0.0, 17M26D:2:0.0, 36M1D:2:0.0, 14M27D:2:0.0, 29M12D:2:0.0, 27M9I:2:0.0, 33M17D:2:0.0, 19M25D:2:0.0, 9M20D:2:0.0, 36M16D:2:0.0, 30M3D:2:0.0, 33M7I:2:0.0, 27M19I:2:0.0, 23M2D:2:0.0, 32M6I:2:0.0, 33M19D:2:0.0, 35M1I:2:0.0, 16M19D:2:0.0, 5M25D:2:0.0, 24M47I:2:0.0, 30M2I:2:0.0, 6M17D:2:0.0, 28M13I:2:0.0, 22M24D:2:0.0, 27M5D:2:0.0, 20M3D:2:0.0, 28M9D:2:0.0, 28M6I:2:0.0, 11M16D:2:0.0, 19M16D:2:0.0, 30M27I:2:0.0, 24M11I:2:0.0, 19M22D:2:0.0, 24M8D:2:0.0, 14M22D:2:0.0, 31M10D:2:0.0, 10M23D:2:0.0, 28M10I:2:0.0, 15M19D:2:0.0, 22M18I:2:0.0, 17M24D:2:0.0, 33M3I:2:0.0, 35M4I:2:0.0, 6M21D:2:0.0, 33M16D:2:0.0, 30M4I:2:0.0, 26M2D:2:0.0, 28M18I:2:0.0, 17M16D:2:0.0, 12M17D:1:0.0, 32M4D:1:0.0, 31M9I:1:0.0, 13M22D:1:0.0, 30M10I:1:0.0, 9M27D:1:0.0, 28M5D:1:0.0, 10M15D:1:0.0, 35M18D:1:0.0, 20M22D:1:0.0, 29M15I:1:0.0, 31M5D:1:0.0, 24M8I:1:0.0, 24M30D:1:0.0, 14M30D:1:0.0, 24M26D:1:0.0, 34M18I:1:0.0, 10M14D:1:0.0, 30M8D:1:0.0, 26M21D:1:0.0, 29M27I:1:0.0, 24M12I:1:0.0, 22M2D:1:0.0, 20M19D:1:0.0, 19M31D:1:0.0, 32M12I:1:0.0, 31M1I:1:0.0, 28M8D:1:0.0, 26M10I:1:0.0, 20M26D:1:0.0, 29M13I:1:0.0, 25M20D:1:0.0, 23M14D:1:0.0, 35M13D:1:0.0, 24M9I:1:0.0, 14M16D:1:0.0, 20M30D:1:0.0, 7M34D:1:0.0, 31M9D:1:0.0, 27M1I:1:0.0, 4M43D:1:0.0, 31M19D:1:0.0, 32M1I:1:0.0, 10M17D:1:0.0, 15M26D:1:0.0, 29M4D:1:0.0, 23M3I:1:0.0, 4M34D:1:0.0, 29M3I:1:0.0, 32M9I:1:0.0, 27M19D:1:0.0, 25M10D:1:0.0, 27M11I:1:0.0, 29M10I:1:0.0, 5M28D:1:0.0, 29M26I:1:0.0, 26M22D:1:0.0, 32M14D:1:0.0, 24M38I:1:0.0, 31M5I:1:0.0, 29M53I:1:0.0, 27M6I:1:0.0, 33M14D:1:0.0, 6M19D:1:0.0, 32M19D:1:0.0, 4M36D:1:0.0, 13M12D:1:0.0, 23M15I:1:0.0, 30M37I:1:0.0, 36M19D:1:0.0, 13M32D:1:0.0, 24M48I:1:0.0, 29M16I:1:0.0, 31M6D:1:0.0, 26M17D:1:0.0, 14M36D:1:0.0, 30M7I:1:0.0, 22M15I:1:0.0, 15M16D:1:0.0, 26M9D:1:0.0, 11M40D:1:0.0, 12M30D:1:0.0, 24M7D:1:0.0, 12M34D:1:0.0, 29M20I:1:0.0, 23M14I:1:0.0, 29M8D:1:0.0, 29M12I:1:0.0, 25M14D:1:0.0, 22M19I:1:0.0, 30M2D:1:0.0, 15M22D:1:0.0, 20M17D:1:0.0, 11M31D:1:0.0, 25M5I:1:0.0, 28M4D:1:0.0, 24M6I:1:0.0, 35M50I:1:0.0, 22M34I:1:0.0, 14M11D:1:0.0, 28M15I:1:0.0, 24M21D:1:0.0, 26M15D:1:0.0, 31M27I:1:0.0, 25M30I:1:0.0, 24M13I:1:0.0, 27M33I:1:0.0, 16M26D:1:0.0, 4M38D:1:0.0, 30M15I:1:0.0, 11M30D:1:0.0, 23M1D:1:0.0, 25M11D:1:0.0, 29M11D:1:0.0, 34M6D:1:0.0, 26M11I:1:0.0, 30M1D:1:0.0, 33M10D:1:0.0, 31M2I:1:0.0, 36M23I:1:0.0, 17M18D:1:0.0, 20M31D:1:0.0, 4M25D:1:0.0, 14M31D:1:0.0, 6M23D:1:0.0, 33M20I:1:0.0, 12M21D:1:0.0, 24M36I:1:0.0, 35M5I:1:0.0, 31M8D:1:0.0, ")
# test = debug.parse_info()
test = debug.process()

print()

# TODO: Polishing; Count the number of files to process => automatically processing all files in the directory
