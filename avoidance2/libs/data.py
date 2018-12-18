

aa2codon =  {'phe' : ['ttt','ttc'] ,
            'leu' : ['tta','ttg','ctt','ctc','cta','ctg'],
            'ile' : ['att','atc','ata'],
            'met' : ['atg'],
            'val' : ['gtt','gtc','gta','gtg'],
            'ser' : ['tct','tcc','tca','tcg','agt','agc'],
            'pro' : ['cct','ccc','cca','ccg'],
            'thr' : ['act','acc','aca','acg'],
            'ala' : ['gct','gcc','gca','gcg'],
            'tyr' : ['tat','tac'],
            'stop': ['taa','tag','tga'],
            'his' : ['cat','cac'],
            'gln' : ['caa','cag'],
            'asn' : ['aat','aac'],
            'lys' : ['aaa','aag'],
            'asp' : ['gat','gac'],
            'glu' : ['gaa','gag'],
            'cys' : ['tgt','tgc'],
            'trp' : ['tgg'],
            'arg' : ['cgt','cgc','cga','cgg','aga','agg'],
            'gly' : ['ggt','ggc','gga','ggg'] }


codon2aa =  { 'aaa' : 'lys',
            'aag' : 'lys',
            'tta' : 'leu',
            'ttg' : 'leu',
            'ctt' : 'leu',
            'ctc' : 'leu',
            'cta' : 'leu',
            'ctg' : 'leu',
            'gaa' : 'glu',
            'gag' : 'glu',
            'caa' : 'gln',
            'cag' : 'gln',
            'cat' : 'his',
            'cac' : 'his',
            'aat' : 'asn',
            'aac' : 'asn',
            'tct' : 'ser',
            'tcc' : 'ser',
            'tca' : 'ser',
            'tcg' : 'ser',
            'agt' : 'ser',
            'agc' : 'ser',
            'cgt' : 'arg',
            'cgc' : 'arg',
            'cga' : 'arg',
            'cgg' : 'arg',
            'aga' : 'arg',
            'agg' : 'arg',
            'tgg' : 'trp',
            'gct' : 'ala',
            'gcc' : 'ala',
            'gca' : 'ala',
            'gcg' : 'ala',
            'tgt' : 'cys',
            'tgc' : 'cys',
            'ggt' : 'gly',
            'ggc' : 'gly',
            'gga' : 'gly',
            'ggg' : 'gly',
            'gat' : 'asp',
            'gac' : 'asp',
            'ttt' : 'phe',
            'ttc' : 'phe',
            'atg' : 'met',
            'tat' : 'tyr',
            'tac' : 'tyr',
            'gtt' : 'val',
            'gtc' : 'val',
            'gta' : 'val',
            'gtg' : 'val',
            'act' : 'thr',
            'acc' : 'thr',
            'aca' : 'thr',
            'acg' : 'thr',
            'cct' : 'pro',
            'ccc' : 'pro',
            'cca' : 'pro',
            'ccg' : 'pro',
            'att' : 'ile',
            'atc' : 'ile',
            'ata' : 'ile',
            'taa' : 'stop',
            'tag' : 'stop',
            'tga' : 'stop' }


codon_to_n = {'aaa': 0,
             'aac': 1,
             'aag': 2,
             'aat': 3,
             'aca': 4,
             'acc': 5,
             'acg': 6,
             'act': 7,
             'aga': 8,
             'agc': 9,
             'agg': 10,
             'agt': 11,
             'ata': 12,
             'atc': 13,
             'atg': 14,
             'att': 15,
             'caa': 16,
             'cac': 17,
             'cag': 18,
             'cat': 19,
             'cca': 20,
             'ccc': 21,
             'ccg': 22,
             'cct': 23,
             'cga': 24,
             'cgc': 25,
             'cgg': 26,
             'cgt': 27,
             'cta': 28,
             'ctc': 29,
             'ctg': 30,
             'ctt': 31,
             'gaa': 32,
             'gac': 33,
             'gag': 34,
             'gat': 35,
             'gca': 36,
             'gcc': 37,
             'gcg': 38,
             'gct': 39,
             'gga': 40,
             'ggc': 41,
             'ggg': 42,
             'ggt': 43,
             'gta': 44,
             'gtc': 45,
             'gtg': 46,
             'gtt': 47,
             'tac': 48,
             'tat': 49,
             'tca': 50,
             'tcc': 51,
             'tcg': 52,
             'tct': 53,
             'tgc': 54,
             'tgg': 55,
             'tgt': 56,
             'tta': 57,
             'ttc': 58,
             'ttg': 59,
             'ttt': 60}

n_to_codon = {0: 'aaa',
            1: 'aac',
            2: 'aag',
            3: 'aat',
            4: 'aca',
            5: 'acc',
            6: 'acg',
            7: 'act',
            8: 'aga',
            9: 'agc',
            10: 'agg',
            11: 'agt',
            12: 'ata',
            13: 'atc',
            14: 'atg',
            15: 'att',
            16: 'caa',
            17: 'cac',
            18: 'cag',
            19: 'cat',
            20: 'cca',
            21: 'ccc',
            22: 'ccg',
            23: 'cct',
            24: 'cga',
            25: 'cgc',
            26: 'cgg',
            27: 'cgt',
            28: 'cta',
            29: 'ctc',
            30: 'ctg',
            31: 'ctt',
            32: 'gaa',
            33: 'gac',
            34: 'gag',
            35: 'gat',
            36: 'gca',
            37: 'gcc',
            38: 'gcg',
            39: 'gct',
            40: 'gga',
            41: 'ggc',
            42: 'ggg',
            43: 'ggt',
            44: 'gta',
            45: 'gtc',
            46: 'gtg',
            47: 'gtt',
            48: 'tac',
            49: 'tat',
            50: 'tca',
            51: 'tcc',
            52: 'tcg',
            53: 'tct',
            54: 'tgc',
            55: 'tgg',
            56: 'tgt',
            57: 'tta',
            58: 'ttc',
            59: 'ttg',
            60: 'ttt'}



if __name__ == '__main__':
    pass

