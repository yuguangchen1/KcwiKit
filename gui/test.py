import tkinter as tk
from tkinter import messagebox

def update_segments():
    try:
        # Get the content of the Text widget
        content = text_widget.get("1.0", tk.END).strip()
        
        # Split the content by lines
        lines = content.split('\n')
        
        # Temporary lists to hold new segments and widths
        new_segments = []
        new_widths = []
        
        for line in lines:
            # Split each line by comma and convert to integers/floats
            start, end, width = line.split(',')
            start = int(start.strip())
            end = int(end.strip())
            width = float(width.strip())
            
            # Append to temporary lists
            new_segments.append((start, end, width))
        
        # Update global segments and widths
        global segments, widths
        segments = [(start, end) for start, end, width in new_segments]
        widths = [width for _, _, width in new_segments]
        
        messagebox.showinfo("Update", "Segments and widths have been updated successfully!")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to update segments and widths: {e}")

# Initial segments and widths
segments = [(0, 1000), (1000, 5000)]
widths = [50, 20]

# Create the main window
root = tk.Tk()
root.title("Edit Segments")

# Create a Text widget
text_widget = tk.Text(root, height=10, width=50)
text_widget.pack()

# Insert initial segment contents into the Text widget
for (start, end), width in zip(segments, widths):
    text_widget.insert(tk.END, f'{start}, {end}, {width}\n')

# Create an Update button
update_button = tk.Button(root, text="Update Segments", command=update_segments)
update_button.pack()

# Start the Tkinter main loop
root.mainloop()
