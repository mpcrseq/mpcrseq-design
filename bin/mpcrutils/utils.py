import json



def parse_mfeprimer(jsoninfile):
    with open(jsoninfile) as f:
      mfepdata = json.load(f)
      