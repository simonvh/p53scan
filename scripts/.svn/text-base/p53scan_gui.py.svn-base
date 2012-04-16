#!/usr/bin/env python
import sys
import string
import os
from Tkinter import *
import tkFileDialog
import tkSimpleDialog
import tkMessageBox
from PftScan.Fasta import Fasta
from PftScan.core import scan

# Constants
NAME = "p53scan"
VERSION = "1.05"
DESCRIPTION = """%s - Version %s
p53scan (http://www.ncmls.nl/bioinfo/p53scan/) is an algorithm to locate
p53 binding sites in DNA sequences. It is described in:
Smeenk L, van Heeringen SJ, Koeppel M, Driel MA, Bartels SJ, Akkers RC, 
Denissov S, Stunnenberg HG, Lohrum M. Characterization of genome-wide 
p53-binding sites upon stress response. Nucleic Acids Res. (2008) """ % (NAME, VERSION)

# Constants determining default values for p53scan
DEFAULT_MIN_SPACER = 0
DEFAULT_MAX_SPACER = 13
DEFAULT_NREPORT = 1 
DEFAULT_CUTOFFS = [4.393,13.7,14.7,15.9,14.8,13.8,15.2,15.6,13.4,14.3,13.0,14.7,13.4,15.5] 

PWM_LEFT = [
[0.3611,0.0769,0.4003,0.1664],
[0.2716,0.0283,0.5667,0.1381],
[0.6358,0.0016,0.3344,0.0330],
[0.0016,0.9859,0.0016,0.0157],
[0.8085,0.0063,0.0502,0.1397],
[0.2276,0.0157,0.0330,0.7284],
[0.0031,0.0016,0.9984,0.0016],
[0.0377,0.3799,0.0016,0.5856],
[0.0816,0.7096,0.0173,0.1962],
[0.1350,0.4035,0.0675,0.3987],
]

PWM_RIGHT = [
[0.4035,0.0534,0.4129,0.1350],
[0.1743,0.0126,0.7394,0.0785],
[0.6044,0.0016,0.3673,0.0314],
[0.0016,0.9984,0.0016,0.0031],
[0.6970,0.0298,0.0188,0.2590],
[0.1381,0.0471,0.0078,0.8116],
[0.0173,0.0031,0.9827,0.0016],
[0.0377,0.3501,0.0031,0.6138],
[0.1413,0.5385,0.0235,0.3014],
[0.1570,0.4066,0.0911,0.3501],
]

class ScanDialog(tkSimpleDialog.Dialog):
	def switch_default(self):
		if self.default_switches[0]['state'] == DISABLED:
			for entry in self.default_switches:
				entry['state'] = NORMAL
			if self.r.get() == 3:
				self.cutoff_entry['state'] = NORMAL
		else:
			for entry in self.default_switches + self.entries:
				entry['state'] = DISABLED

	def cutoff_default(self):
		self.cutoff_entry['state'] = NORMAL
		self.cutoff_entry.delete(0,END)
		self.cutoff_entry.insert(END, "Default")
		self.cutoff_entry['state'] = DISABLED

	def cutoff_none(self):
		self.cutoff_entry['state'] = NORMAL
		self.cutoff_entry.delete(0,END)
		self.cutoff_entry.insert(END, "None")
		self.cutoff_entry['state'] = DISABLED

	def cutoff_on(self):
		self.cutoff_entry['state'] = NORMAL
		self.cutoff_entry.delete(0,END)
		self.cutoff_entry.insert(END, "0")
		
	def body(self, master):
		self.scan_ok = False

		boxes = [
			("Minspacer:", 0, DEFAULT_MIN_SPACER, True),
			("Maxspacer:", 1, DEFAULT_MAX_SPACER, True),
			("Numreport:", 2, DEFAULT_NREPORT, True),
			("Score cutoff:", 4, "Default", False),
		]
	
		self.entries = []
		self.default_switches = []

		for (label, r, val, ds) in boxes:
			Label(master, text=label).grid(row=r)
			entry = Entry(master)
			entry.grid(row=r, column=1, columnspan=3, sticky=W)
			entry.insert(END, val)
			self.entries.append(entry)
			if ds:
				self.default_switches.append(entry)
		
		self.cutoff_entry = self.entries[-1]

		Label(master).grid(row=3)
		Label(master).grid(row=6)

		self.r = IntVar()
		self.r.set(1)
		for (label, v, col, function) in [("Default", 1, 0, self.cutoff_default), ("None", 2, 1, self.cutoff_none),("Cutoff", 3, 2, self.cutoff_on)]:
			r = Radiobutton(master, text=label, variable=self.r, value=v, command=function)
			r.grid(row=5, column=col, sticky=W)
			self.default_switches.append(r)

		for entry in self.entries + self.default_switches:
			entry["state"] = DISABLED

		self.cb = Checkbutton(master, text="Default Settings", command=self.switch_default)
		self.cb.grid(row=7, columnspan=4, sticky=W)
		self.cb.select()
		
		return self.cb # initial focus

	def apply(self):
		pass

	def validate(self):
		try:
			self.scan_ok = True
			self.min_spacer = string.atoi(self.entries[0].get())
			self.max_spacer = string.atoi(self.entries[1].get())
			self.nreport = string.atoi(self.entries[2].get())
			self.defaults = True
			self.cutoff = None
			if self.r.get() > 1:
				self.defaults = False
				if self.r.get() == 3:
					self.cutoff = string.atof(self.entries[3].get())
			return 1
		except:
			tkMessageBox.showwarning("Bad input", "Illegal values, please try again")
			return 0

class MultiListbox(Frame):
	# MultiListbox Tkinter widget code by Brent Burley
	# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52266
	# SortableTable widget additions by Rick Lawson
	# http://tkinter.unpythonic.net/wiki/SortableTable
	def __init__(self, master, lists):
		Frame.__init__(self, master)
		self.lists = []
		self.colmapping={}
		self.origData = None
		for l,w in lists:
			frame = Frame(self); frame.pack(side=LEFT, expand=YES, fill=BOTH)
			b = Button(frame, text=l, borderwidth=1, relief=RAISED)
			b.pack(fill=X)
			b.bind('<Button-1>', self._sort)
			self.colmapping[b]=(len(self.lists),1)
			lb = Listbox(frame, width=w, borderwidth=0, selectborderwidth=0,
						 relief=FLAT, exportselection=FALSE)
			lb.pack(expand=YES, fill=BOTH)
			self.lists.append(lb)
			lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))
			lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
			lb.bind('<Leave>', lambda e: 'break')
			lb.bind('<B2-Motion>', lambda e, s=self: s._b2motion(e.x, e.y))
			lb.bind('<Button-2>', lambda e, s=self: s._button2(e.x, e.y))
		frame = Frame(self); frame.pack(side=LEFT, fill=Y)
		Label(frame, borderwidth=1, relief=RAISED).pack(fill=X)
		sb = Scrollbar(frame, orient=VERTICAL, command=self._scroll)
		sb.pack(expand=YES, fill=Y)
		self.lists[0]['yscrollcommand']=sb.set

	def _select(self, y):
		row = self.lists[0].nearest(y)
		self.selection_clear(0, END)
		self.selection_set(row)
		return 'break'

	def _button2(self, x, y):
		for l in self.lists: l.scan_mark(x, y)
		return 'break'

	def _b2motion(self, x, y):
		for l in self.lists: l.scan_dragto(x, y)
		return 'break'

	def _scroll(self, *args):
		for l in self.lists:
			apply(l.yview, args)

	def curselection(self):
		return self.lists[0].curselection()

	def delete(self, first, last=None):
		for l in self.lists:
			l.delete(first, last)

	def get(self, first, last=None):
		result = []
		for l in self.lists:
			result.append(l.get(first,last))
		if last:
			if last - first == 0:
				return [result]
			else:
				return apply(map, [None] + result)
		return [result]

	def index(self, index):
		self.lists[0].index(index)

	def insert(self, index, *elements):
		for e in elements:
			i = 0
			for l in self.lists:
				l.insert(index, e[i])
				i = i + 1

	def size(self):
		return self.lists[0].size()

	def see(self, index):
		for l in self.lists:
			l.see(index)

	def selection_anchor(self, index):
		for l in self.lists:
			l.selection_anchor(index)

	def selection_clear(self, first, last=None):
		for l in self.lists:
			l.selection_clear(first, last)

	def selection_includes(self, index):
		return self.lists[0].selection_includes(index)

	def selection_set(self, first, last=None):
		for l in self.lists:
			l.selection_set(first, last)

	def _sort(self, e):
		# get the listbox to sort by (mapped by the header button)
		b=e.widget
		col, direction = self.colmapping[b]

		# get the entire table data into mem
		tableData = self.get(0,END)
		if self.origData == None:
			import copy
			self.origData = copy.deepcopy(tableData)

		rowcount = len(tableData)

		#remove old sort indicators if it exists
		for btn in self.colmapping.keys():
			lab = btn.cget('text')
			if lab[0]=='[': btn.config(text=lab[4:])

		btnLabel = b.cget('text')
		#sort data based on direction
		if direction==0:
			tableData = self.origData
		else:
			if direction==1: b.config(text='[+] ' + btnLabel)
			else: b.config(text='[-] ' + btnLabel)
			# sort by col
			def colsort(x, y, mycol=col, direction=direction):
				return direction*cmp(x[mycol], y[mycol])

			tableData.sort(colsort)

		#clear widget
		self.delete(0,END)

		# refill widget
		for row in range(rowcount):
			self.insert(END, tableData[row])

		# toggle direction flag 
		if(direction==1): direction=-1
		else: direction += 1
		self.colmapping[b] = (col, direction)
	


class StatusBar(Frame):

	def __init__(self, master):
		Frame.__init__(self, master)
		self.label = Label(self, bd=1, relief=SUNKEN, anchor=W)
		self.label.pack(fill=X)

	def set(self, format, *args):
		self.label.config(text=format % args)
		self.label.update_idletasks()

	def clear(self):
		self.label.config(text="")
		self.label.update_idletasks()

class App:

	def open_fasta_file(self):
		name = tkFileDialog.askopenfilename()
		
		if name:
			try:
				self.f = Fasta(name)
				self.loaded_file = os.path.split(name)[-1]
				self.mlb.delete(0, self.mlb.size() - 1)
				self.status.set("Loaded %s - %s" % (self.loaded_file, self.f))
				self.scanmenu.entryconfig(0, state=NORMAL)
			except:
				
				tkMessageBox.showerror("Error opening file", "Error opening %s.\nIs it a valid FASTA file?" % os.path.split(name)[-1])


	def save_results(self):
		name = tkFileDialog.asksaveasfilename(filetypes=[("gff", "*.gff")], defaultextension=".gff")
		if name:
			f = open(name, "w")	
			for row in self.mlb.get(0, self.mlb.size() - 1):
				f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s %s %s\n" % (row[0], "p53scan", "p53bs", row[1], row[2], row[3], "+", ".", row[4], row[5], row[6]))
			f.close()
			self.status.set("p53scan results written to %s" % os.path.split(name)[-1])

	def scan_file(self):
		d = ScanDialog(self.master)
		if d.scan_ok:
			min_spacer = d.min_spacer
			max_spacer = d.max_spacer
			nreport = d.nreport
			cutoffs = []
			if d.defaults:
				cutoffs = DEFAULT_CUTOFFS
			elif d.cutoff:
				cutoffs = (max_spacer - min_spacer + 1) * [d.cutoff]
			
			d.destroy()

			self.status.set("Scanning %s for p53 binding sites" % self.loaded_file)
			self.mlb.delete(0, self.mlb.size() - 1)
			count = 0
			for (id,seq) in self.f.items():
				result = scan(seq.upper(), PWM_LEFT, PWM_RIGHT, min_spacer, max_spacer, cutoffs, nreport)
				for (score, pos, spacer, strand) in result:
					count += 1
					self.mlb.insert(END,(id, pos, pos + len(PWM_LEFT) + spacer + len(PWM_RIGHT), score,
						seq[pos: pos + len(PWM_LEFT)],
						seq[pos + len(PWM_LEFT): pos + len(PWM_LEFT) + spacer],
						seq[pos + len(PWM_LEFT) + spacer: pos + len(PWM_LEFT) + spacer + len(PWM_RIGHT)]
					))
			self.status.set("Done scanning %s. Found %s binding sites." % (self.loaded_file, count))
			self.filemenu.entryconfig(1, state=NORMAL)

	def exit(self):
		sys.exit()

	def resize(self, event):
		frame.pack()

        def about(self):
            tkMessageBox.showinfo("About %s" % NAME, DESCRIPTION)

	def __init__(self, master): 
		# create a menu
                self.master = master

                master.title("p53scan")
		menu = Menu(master)
		master.config(menu=menu)

		self.filemenu = Menu(menu,tearoff=0)
		menu.add_cascade(label="File", menu=self.filemenu)
		self.filemenu.add_command(label="Open Fastafile...", command=self.open_fasta_file)
		self.filemenu.add_command(label="Save Results...", command=self.save_results)
		self.filemenu.add_separator()
		self.filemenu.add_command(label="Exit", command=self.exit)
                self.filemenu.entryconfig(1, state=DISABLED)

		self.scanmenu = Menu(menu,tearoff=0)
		menu.add_cascade(label="Scan", menu=self.scanmenu)
		self.scanmenu.add_command(label="Scan file", command=self.scan_file)
                self.scanmenu.entryconfig(1, state=DISABLED)

		helpmenu = Menu(menu,tearoff=0)
		menu.add_cascade(label="Help", menu=helpmenu)
		helpmenu.add_command(label="About...", command=self.about)

		self.status = StatusBar(master)
		self.status.pack(side=BOTTOM, fill=X)
		self.status.set("No file loaded")

		self.mlb = MultiListbox(master, (('id', 30), ('start', 4),('end',4), ('score', 4), ('left', 15), ('spacer', 10), ('right', 15)))
		self.mlb.pack(expand=YES, fill=BOTH)

root = Tk()
app = App(root)

mainloop()

