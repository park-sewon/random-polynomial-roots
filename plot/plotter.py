import matplotlib.pyplot as plt
import matplotlib.patches as patches
def add_by_path(ax1, s, w):
    area_x = [-w/2,w/2]
    area_y = [-w/2,w/2]
    # print s
    for o in range(1,len(s)):
        if s[o] == "1":
            area_x[0] = (area_x[0] + area_x[1]) / 2
            area_y[0] = (area_y[0] + area_y[1]) / 2
        elif s[o] == "2":
            area_x[0] = (area_x[0] + area_x[1]) / 2
            area_y[1] = (area_y[0] + area_y[1]) / 2
        elif s[o] == "3":
            area_x[1] = (area_x[0] + area_x[1]) / 2
            area_y[1] = (area_y[0] + area_y[1]) / 2
        elif s[o] == "4":
            area_x[1] = (area_x[0] + area_x[1]) / 2
            area_y[0] = (area_y[0] + area_y[1]) / 2

    # print area_x
    # print area_y
    ax1.add_patch(
        patches.Rectangle(
            (area_x[0], area_y[0]),  # (x,y)
            area_x[1] - area_x[0],  # width
            area_y[1] - area_y[0],  # height
        fill=False, edgecolor='b')
    )
    return ax1

def add_by_frac(ax1, f):
    c_re = float(f.split(",")[0])
    c_im = float(f.split(",")[1])
    w = float(f.split(",")[2])
    ax1.add_patch(
        patches.Rectangle(
            (c_re - w/2, c_im-w/2),  # (x,y)
            w,  # width
            w,  # height
        fill=False, edgecolor='b')
    )
    return ax1

if __name__=="__main__":
    with open("data/newton2.txt","r") as f:



        fig1 = plt.figure()
        initial_width = 4.0
        fig1.set_figheight(30)
        fig1.set_figwidth(30)
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.xlim([-initial_width/2, initial_width/2])
        plt.ylim([-initial_width/2, initial_width/2])
        for line in f:
            line = line.replace("\n","")
            ax1 = add_by_frac(ax1, line)
    fig1.savefig('img/newtdon'+str(10)+'.png', dpi=600, bbox_inches='tight')


    s = "-.556274414063E-0003"
    f = float(s)
    print f
    quit(0)
    initial_width = 8.0

    dic = {}
    with open("data/gauss_10.txt","r") as f:
        for line in f:
            line = line.replace("\n","")
            if dic.has_key(len(line)):
                dic[len(line)].append(line)
            else:
                dic[len(line)] = [line]
    fig1 = plt.figure()

    for v in range(1,len(dic)):
        print dic[v+1]
        fig1.set_figheight(15)
        fig1.set_figwidth(15)
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.xlim([-initial_width/2, initial_width/2])
        plt.ylim([-initial_width/2, initial_width/2])
        for s in dic[v+1]:
            ax1 = add_by_path(ax1, s,initial_width)
    fig1.savefig('img/newton'+str(v)+'.png', dpi=300, bbox_inches='tight')
