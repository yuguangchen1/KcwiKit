import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

class InteractivePlot(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Interactive Plot")
        self.geometry("800x600")

        self.frame = tk.Frame(self)
        self.frame.pack(fill=tk.BOTH, expand=True)

        # Matplotlib figure
        self.figure = plt.figure()
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Example data
        self.data = [(10, 50), (50, 80), (90, 30), (130, 70), (170, 40),
                     (210, 90), (250, 60), (290, 40), (330, 100), (370, 80)]

        self.selected_regions = []
        self.mode = None
        self.start_point = None

        self.bind("<KeyPress-e>", lambda event: self.set_mode("exclude"))
        self.bind("<KeyPress-I>", lambda event: self.set_mode("include"))

        self.canvas.mpl_connect('button_press_event', self.on_click)

        self.draw_data()

    def set_mode(self, mode):
        self.mode = mode

    def on_click(self, event):
        if self.mode == "exclude" or self.mode == "include":
            if self.start_point is None:
                self.start_point = (event.xdata, event.ydata)
            else:
                end_point = (event.xdata, event.ydata)
                self.selected_regions.append((self.start_point, end_point))
                self.draw_region(self.start_point, end_point, self.mode)
                self.start_point = None

    def draw_data(self):
        x_data, y_data = zip(*self.data)
        self.ax.plot(x_data, y_data, 'bo-')
        self.canvas.draw()

    def draw_region(self, start_point, end_point, mode):
        if mode == "exclude":
            color = "red"
        else:
            color = "green"
        rect = plt.Rectangle(start_point, end_point[0] - start_point[0], end_point[1] - start_point[1],
                             edgecolor=color, facecolor=color, alpha=0.3)
        self.ax.add_patch(rect)
        self.canvas.draw()

if __name__ == "__main__":
    app = InteractivePlot()
    app.mainloop()
