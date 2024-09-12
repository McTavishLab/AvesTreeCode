import csv
"""
A function to translate crosswalk CSV tables into python dictionary objects.
returns a dictionary with ott_ids as keys and 'alt_name' as values.
"""
def crosswalk_to_dict(taxon_crosswalk_file,
                      ott_id_key="ott_id",
                      alt_name='SCI_NAME',
                      tax_filter='species', 
                      delimiter=","):
                      name_map = {}
                      no_maps = []
                      double_maps = {}
                      with open(taxon_crosswalk_file, newline='') as csvfile:
                        taxreader = csv.reader(csvfile, delimiter=delimiter)
                        i=0
                        for lii in taxreader:
                            i+=1
                            if i ==1:
                                header = lii
                                ott_id_index = header.index(ott_id_key)
                                alt_name_index = header.index(alt_name)
                            else:
                                ott_id = 'ott'+lii[ott_id_index]
                                if ott_id in name_map:
                                    if ott_id in double_maps:
                                        double_maps[ott_id].append(dict(zip(header, lii)))
                                    else:
                                        double_maps[ott_id ]=[dict(zip(header, lii))]
                                name_map[ott_id] = lii[alt_name_index]
                      return name_map