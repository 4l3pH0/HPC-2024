import matplotlib.pyplot as plt
import sys

y = [1, 0.650, 1.812, 3.284, 5.473]
x = [1, 2, 4, 8, 16]

y2 = [1, 3.684052, 6.575054]
x2 = [1, 4, 8]

y3 = [1, 1.938, 3.754]
x3 = [1, 2, 4]

y4 = [1, 2.635, 3.859]
x4 = [1, 2, 4]

y5 = [1, 0.587, 1.682, 2.932, 3.386]
x5 = [1, 2, 4, 8, 16]

y6 = [1, 1.730, 3.431, 2.714]
x6 = [1, 4, 8, 16]


if sys.argv[1] == "80":
    plt.plot(x,y)

    plt.ylabel('Ускорение')
    plt.xlabel('Число нитей')
    plt.title("Ускорение на сетке 80x90")
    plt.show()
elif sys.argv[1] == "160":
    plt.plot(x2,y2)

    plt.ylabel('Ускорение')
    plt.xlabel('Число нитей')
    plt.title("Ускорение на сетке 160x180")
    plt.show()
elif sys.argv[1] == "80m":
    plt.plot(x3,y3)
    plt.ylabel('Ускорение')
    plt.xlabel('Число процессов')
    plt.title("Ускорение на сетке 80x90")
    plt.show()
elif sys.argv[1] == "160m":
    plt.plot(x4,y4)
    plt.ylabel('Ускорение')
    plt.xlabel('Число процессов')
    plt.title("Ускорение на сетке 160x180")
    plt.show()
elif sys.argv[1] == "90mp":
    plt.plot(x5,y5)
    plt.ylabel('Ускорение')
    plt.xlabel('Число процессов x число нитей')
    plt.title("Ускорение на сетке 80x90")
    plt.show()
elif sys.argv[1] == "160mp":
    plt.plot(x6,y6)
    plt.ylabel('Ускорение')
    plt.xlabel('Число процессов x число нитей')
    plt.title("Ускорение на сетке 160x180")
    plt.show()
