source(here::here('data', 'nork_data_2022.r'))
# inputs ------------------
data <- list(ages = ages,
             years = years,
             length_bins = 15:45,
             waa = waa,
             wt_mature = maa * waa / 2,
             spawn_mo = 5,
             catch_obs = catch_obs,
             catch_wt = c(rep(5, 17), rep(50, 45)),
             srv_obs = srv_obs,
             srv_ind = srv_ind,
             srv_sd = srv_sd,
             srv_wt = 0.25,
             fish_age_obs = fish_age_obs,
             fish_age_ind = fish_age_ind,
             fish_age_iss = fish_age_iss,
             fish_age_wt = 0.5,
             srv_age_obs = srv_age_obs,
             srv_age_ind = srv_age_ind,
             srv_age_iss = srv_age_iss,
             srv_age_wt = 0.5,
             fish_size_obs = fish_size_obs,
             fish_size_ind = fish_size_ind,
             fish_size_iss = fish_size_iss,
             fish_size_wt = 0.5,
             age_error = age_error,
             size_age = size_age,
             wt_fmort_reg = 0.1,
             mean_M = 0.06,
             sd_M = 0.05,
             mean_q = 1,
             sd_q = 0.45,
             mean_sigmaR = 1.5,
             sd_sigmaR = 0.01,
             proj_rec_yrs = proj_rec_yrs,
             yield_ratio = yield_ratio)
# parameters ----------------
pars <- list(log_M = log(0.05949983),
             log_a50C = log(8.23719578492),
             deltaC = 1.91866640868,
             log_a50S = log(9.09355629491),
             deltaS = 4.31923047888,
             log_q = log(0.8649374),
             log_mean_R = 3.50391,
             init_log_Rt = c(-1.31784331510, -1.41353004574, -1.45443795306, -1.44801600013, -1.44317223143,
                             -1.47227774217, -1.51491209802, -1.54201844034, -1.54900822105, -1.54357308125, -1.53271591530,
                             -1.51992650551, -1.50655074209, -1.49298490737, -1.47936419813, -1.46579792964, -1.45236085259,
                             -1.43911454732, -1.42609004023, -1.41334175243, -1.40089855851, -1.38878689192, -1.37702228101,
                             -1.36561634960, -1.35457991448, -1.34391631981, -1.33362760589, -1.32371409445, -1.31417494445,
                             -1.30500561365, -1.29620121156, -1.28775584994, -1.27966176612, -1.27190998138, -1.26449234040,
                             -1.25740132870, -1.25062646041, -1.24415861470, -1.23798650648, -1.23210097519, -1.22649216166,
                             -1.22115043981, -1.21606506031, -1.21122735382, -1.20662780636, -1.20225344500, -1.19809591673,
                             -1.19414750188),
             log_Rt = c(-1.18621056605,
                        -1.03273909462,-0.999797306673 ,-1.17237515921,-1.31481440189, -1.29947684085,
                        -1.19650952919, -1.07269785625, -0.858740244220 ,-0.648271337659, -0.687069961675,
                        0.445219166778, -0.468276695996, -0.563352814004, -0.453027729654, -0.584195538322,
                        -0.429755690388, 0.614844946205, 0.243664403047, -0.400119220417, -0.500707835112,
                        -0.335760494094, -0.00147160235074, 0.205169356072, -0.589731491602, 0.662087934070,
                        -0.0974724220631, -0.772529863247, -0.459399034102, -0.338932984917 ,-1.22764801352,
                        -0.470463642208, -0.821606342481, -0.909630584967, -1.23654078531, 0.436194651298,
                        0.0259279936320, -0.463748734096, -0.359703385166, 0.182589967526, -0.706622061830,
                        -1.01841380990, -0.793317388193, -1.45663999946, -2.24234182204, -2.25003316068,
                        -1.76689920898, -1.93202614806, -1.74725802806, -1.63462815593, -1.62908024719,
                        -1.84514613721, -2.16498687286, -1.76579171915, -1.89444424903, -1.60259969329,
                        -1.75681180416, -1.48983850684, -1.39173905744, -1.25892345656, -1.18461587351, -1.12500000434),
             log_mean_F = -3.58385236595,
             log_Ft =  c(-1.28306124683, 0.140199833981, 0.952641543234, 1.67793968633, 2.28974699704, 1.91823978977,
                         1.51557034959, 1.40878774713, 1.12294109172, 0.646529313912, 1.20801861960, 1.24452195275, 0.942561137836,
                         0.835801151139, 0.798114766724, 0.644693350645, -0.735539824209, -1.02062427087, -0.991743165230,
                         -0.920367489050, -0.403644808221, 0.516721388258, 0.389882559033, -0.981570055502, -2.79223891420,
                         -2.60309212227, -2.00728646156, -1.23024028545, -0.956209854061, -0.890429678442, 0.0349559024681,
                         0.546156698362, 0.0400591009128, 0.234941837688, 0.181355511206, -0.335989512042, -0.463551439451,
                         -0.416583500320, 0.167372898873, -0.293055078553, -0.341087697776, -0.279675928688, 0.155952440289,
                         0.0591776963077, -0.00563314517619, 0.0822175014391, -0.0885440790171, -0.113680406335, -0.120666350781,
                         -0.106971346707, -0.196222087322, 0.238469258157, 0.258217566983, 0.185811847296, 0.161828748432,
                         0.0782651773294, -0.500389551457, -0.207598679635, -0.00551620573419, -0.0956412655191, -0.049808327, -0.241028829),
             log_F35 = log(0.0736156723516),
             log_F40 = log(0.0612999005103),
             log_F50 = log(0.0428699169247),
             sigmaR = 1.5)

# map ----
map = list(log_M = factor(NA),
           log_a50C = factor(NA),
           deltaC = factor(NA),
           log_a50S = factor(NA),
           deltaS = factor(NA),
           log_q = factor(NA),
           log_mean_R = factor(NA),
           init_log_Rt = factor(rep(NA, 48)),
           log_Rt = factor(rep(NA, 62)),
           log_mean_F = factor(NA),
           log_Ft = factor(rep(NA, 62)),
           log_F35 = factor(NA),
           log_F40 = factor(NA),
           log_F50 = factor(NA),
           sigmaR = factor(NA))

out = list(data, pars, map)
