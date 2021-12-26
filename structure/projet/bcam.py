from modeller import *
from modeller.automodel import *


env = Environ()
a = AutoModel(env, alnfile='align.ali',
              knowns='6JXR', sequence='BCAM',
              assess_methods=(assess.DOPE,
                              assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()
