from options import CamSolveOptions, calculate_cam_solve

if __name__ == '__main__':
    cso = CamSolveOptions(kulachok_type = 'thin',
                          N = 1000,
                          graphs_kulachok_flag = True,
                          graphs_profil_flag = True)
    calculate_cam_solve(cso)
