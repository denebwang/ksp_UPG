
from colorama import init

from UPG import Status

init()

####### helper functions ########
def move_cursor(x, y): 
    print("\x1b[{};{}H".format(y, x), end="")
    
def distance_unit_convert(dis):
    if dis > 1e7: 
        return "{:,.2f} Mm".format(dis / 1e6)
    if dis > 1e5:
        return "{:,.2f} km".format(dis / 1e3)

    return "{:,.2f} m".format(dis)


####### printer functions ########
def background():
    print("┌" + "".ljust(60,'─') + "┐")
    print("│" + "UPG Powered Descent".center(60) + "│")
    print("├" + "".ljust(60,'─') + "┤")
    print("│" + "Status: ".rjust(26).ljust(60) + "│")
    print("├" + "".ljust(33,'─') + "┬" + "".ljust(26,'─') + "┤")
    
    print("│" + "Target downrange: ".ljust(33) + "│" + "t_go: ".ljust(26) + "│")
    print("│" + "Guidance downrange: ".ljust(33) + "│" + "v_go: ".ljust(26) + "│")
    
    print("├" + "".ljust(28,'─') + "┬" + "".ljust(4,'─') + "┴" + "".ljust(26,'─') + "┤")
    
    print("│" + "Convergence: ".ljust(28) + "│" + "Norm: ".ljust(31) + "│")
    print("│" + "Last convergence: ".ljust(28) + "│" + "Last err: ".ljust(31) + "│")
    
    print("├" + "".ljust(28,'─') + "┴" + "".ljust(31,'─') + "┤")
    print("│" + "Thrust err:".ljust(23) + "AP err: ".ljust(15) + "Land err:".ljust(22) + "│")
    print("└" +"".ljust(60,'─') + '┘')


def update(status, tdr, gdr, t_go, t_1_go, t_2_go, v_go, conver, norm, last_err, last_conver, thrust_err, ap_err, land_err): 
    print_status(status)
    print_downrange(tdr, gdr)
    print_time(t_go, t_1_go, t_2_go)
    print_vgo(v_go)
    print_solver(conver, norm, last_err, last_conver)
    print_error(thrust_err, ap_err, land_err)


def update_upg(status, tdr, gdr, t_go, t_1_go, t_2_go, v_go, conver, norm, last_err, last_conver, ap_err, land_err):
    update(status, tdr, gdr, t_go, t_1_go, t_2_go, v_go, conver, norm, last_err, last_conver, None, ap_err, land_err)


def update_terminal(dr, t_go, thrust_err, ap_err):
    update(Status.Terminal, dr, None, t_go, None, None, None, None, None, None, None, thrust_err, ap_err, None)


def print_status(status): 
    move_cursor(28, 4)
    str = "".ljust(20)
    if status == Status.PDI:
        str = "Waiting to descent".ljust(20)
    elif status == Status.Powered_descent:
        str = "Main Guidance".ljust(20)
    elif status == Status.Terminal:
        str = "Terminal Guidance".ljust(20)
    elif status == Status.Finished:
        str = "Finished".ljust(20)
    print(str, end='')


def print_downrange(target_dr, guidance_dr):
    # scale to correct unit
    if target_dr is not None:
        tdr_str = distance_unit_convert(target_dr).rjust(14)
    else : tdr_str = "".rjust(14)
    if guidance_dr is not None:
        gdr_str = distance_unit_convert(guidance_dr).rjust(13)
    else : gdr_str = "".rjust(13)
    move_cursor(20, 6)
    print(tdr_str, end='')
    move_cursor(21, 7)
    print(gdr_str, end='')


def print_time(tgo, t1go, t2go):
    move_cursor(42, 6)
    print("{:.1f} s" .format(tgo).rjust(8), end='')
    move_cursor(51, 6)
    if t1go is not None and t2go is not None:
        if t1go > 0.:
            print("(t1 " + "{:.1f})".format(t1go).rjust(7), end='')
            return
        if t2go > 0.:
            print("(t2 " + "{:.1f})".format(t2go).rjust(7), end='')
            return
    print("".rjust(11), end='')


def print_vgo(v):
    move_cursor(42, 7)
    if v is None:
        print("".ljust(20))
    else:
        print("{:,.2f} m/s".format(v).ljust(20), end='')


def print_solver(conver, norm, last_err, last_conver):
    move_cursor(15, 9)
    if conver is None:
        print("".rjust(15), end='')
    else:
        print(str(conver).rjust(15), end='')
    
    move_cursor(39, 9)
    if norm is None:
        print("".rjust(23), end='')
    else:
        print("{:.4e}".format(norm).rjust(23), end='')
    
    move_cursor(20, 10)
    if last_conver is None:
        print("".rjust(10), end='')
    else:
        print("{:.0f} s ago".format(last_conver).rjust(10), end='')
    
    move_cursor(42, 10)
    if last_err is None:
        print("".rjust(20), end='')
    else:
        if last_err == 1:
            s = "" # don't show converge
        if last_err == 2:
            s = "max iter"
        if last_err == 3:
            s = "xtol too small"
        if last_err == 4:
            s = "bad jac"
        if last_err == 5:
            s = "bad improvement"
        print(s.rjust(20), end='')


def print_error(t, ap, l):
    move_cursor(14, 12)
    if t is None:
        print("".rjust(10), end='')
    else:
        print("{:,.1f} kN".format(t / 1e3).rjust(10), end='')
    move_cursor(32, 12)
    print("{:.1f}°".format(ap).rjust(7), end='')
    move_cursor(51, 12)
    if l is None:
        print("".rjust(11), end='')
    else:
        print("{:,.0f} m".format(l).rjust(11), end='')
    
def end():
    move_cursor(1,14)