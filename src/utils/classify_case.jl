#TODO - do we want to change anything here now?
is_baro_wave(params) = all((
    params["config"] == "sphere",
    params["forcing"] == nothing,
    params["surface_scheme"] == nothing,
    params["perturb_initstate"] == true,
))

is_solid_body(params) = all((
    params["config"] == "sphere",
    params["forcing"] == nothing,
    params["rad"] == nothing,
    params["perturb_initstate"] == false,
))

is_column_without_edmf(params) = all((
    params["config"] == "column",
    params["turbconv"] == nothing,
    params["forcing"] == nothing,
    params["turbconv"] != "edmf",
))

is_column_edmf(params) = all((
    params["config"] == "column",
    params["energy_name"] == "rhoe",
    params["forcing"] == nothing,
    params["turbconv"] == "edmf",
    params["rad"] == "DYCOMS_RF01" ||
    params["rad"] == "TRMM_LBA" ||
    params["rad"] == nothing,
))
