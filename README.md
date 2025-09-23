# toolkit

## baraffe_isochrone.py
ages, masses, data_dict = parse_track_file("COND03_models.txt")
interps = build_interpolators(ages, masses, data_dict)
#
Teff_val = interps["Teff"](([5.1, 5.3]))  # 5.1 Myr, 5.3 Mjup
logg_val = interps["g"](([5.1, 5.3]))
print(Teff_val,logg_val)
