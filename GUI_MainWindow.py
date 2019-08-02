from tkinter import *
import ProjectUtils as PJU
import tkinter as tk
import random 
#import pyperclip

class Interface(Frame):
    
    """Main Window.
    All widgets all stocked as attributs in this window."""
    
    def __init__(self, fenetre, **kwargs):
        Frame.__init__(self, fenetre, width=100, height=200, **kwargs)
        self.pack(fill=BOTH)
        self.file_index = 0 #For name file indexing
        #Create widgets
        #
        self.label_windowName = Label(self, text='Hello !', fg='red')
        self.label_windowName.place(relx=0.5, rely=0.0)
        
        self.button_help = Button(self, text='?', bg='blue', fg='white',command=self.doc_Functions)
        self.button_help.place(relx=1.0, rely=0.0, anchor=NE)
        self.help_bool = False
        
        self.cadre_result = LabelFrame(text='result ', width=300, height=500)
        self.cadre_result.place(relx=0.668, rely=0.5, anchor=CENTER)
        
        self.frame_randomDNASeq = Frame(fenetre, width=20, height=10)
        self.frame_randomDNASeq.place(relx=0.16, rely=0.17, anchor=CENTER) 
        self.entryText = StringVar() #To initialize the field to ''
        self.entry_randomDNASeq = Entry(self.frame_randomDNASeq, width=5, textvariable=self.entryText)
        self.entry_randomDNASeq.pack(side='right')
        
        self.button_randomDNASequence = Button(self.frame_randomDNASeq, text='Create a random sequence', command=lambda:self.doc_or_operation(1))
        self.button_randomDNASequence.pack(side='left')
        
        self.button_countAcidNucl = Button(text='Count Nucleic Acid', command=lambda:self.doc_or_operation(3))
        self.button_countAcidNucl.place(relx=0.1, rely=0.24, anchor=CENTER)
           
        self.button_sequenceIsValid = Button(text='Sequence validation', command=lambda:self.doc_or_operation(2))
        self.button_sequenceIsValid.place(relx=0.105, rely=0.31, anchor=CENTER)
        
        self.button_Translate_DNAtoRNA = Button(text='Translate DNA to RNA', command=lambda:self.doc_or_operation(4))
        self.button_Translate_DNAtoRNA.place(relx=0.111, rely=0.380, anchor=CENTER)
        
        self.button_Translate_toProtein = Button(text='Transcribe DNA/RNA to Protein', command=lambda:self.doc_or_operation(5))
        self.button_Translate_toProtein.place(relx=0.145, rely=0.45, anchor=CENTER)
        
        self.button_ReverseDNAComplement = Button(text='Get the DNA reverse Complement', command=lambda:self.doc_or_operation(6))
        self.button_ReverseDNAComplement.place(relx=0.159, rely=0.520, anchor=CENTER)
        
        self.button_Rate_GC = Button(text='Calculate GC rate in sequence', command=lambda:self.doc_or_operation(7))
        self.button_Rate_GC.place(relx=0.144, rely=0.59, anchor=CENTER)
 
        self.button_codonsFrequence = Button(text='Calculate codons frequences of DNA/RNA sequence', command=lambda:self.doc_or_operation(8))
        self.button_codonsFrequence.place(relx=0.225, rely=0.66, anchor=CENTER)
        
        self.button_monoisotopicMass = Button(text='Calculate monoisotopic Mass of sequence', command=lambda:self.doc_or_operation(9))
        self.button_monoisotopicMass.place(relx=0.19, rely=0.73, anchor=CENTER)
  
        self.button_RNASplicing = Button(text='Splice RNA sequence', command=lambda:self.doc_or_operation(10))
        self.button_RNASplicing.place(relx=0.111, rely=0.80, anchor=CENTER)
    
        self.button_DNAAssembler = Button(text='Launch DNA Assembler', command=lambda:self.doc_or_operation(11))
        self.button_DNAAssembler.place(relx=0.12, rely=0.87, anchor=CENTER)
    
        self.button_quit = Button(text='Exit', command=self.quit)
        self.button_quit.place(relx=1.0, rely=1.0, anchor=SE)
        
        self.frame_result = Frame(self.cadre_result, width=500, height=800)   
        self.frame_result.pack(side=RIGHT)
        
        self.scrollbar_cadreResult = Scrollbar(self.frame_result)
        self.scrollbar_cadreResult.pack(side=RIGHT, fill=Y)
        
        self.result_text = Text(self.frame_result, height=30, width=45)
        self.result_text.pack(side=LEFT, fill=Y)
        self.result_text.config(state='disable')
        
        self.scrollbar_cadreResult.config(command=self.result_text.yview)
        self.result_text.config(yscrollcommand=self.scrollbar_cadreResult.set)
        
        self.button_save = Button(text='Save to file', command=lambda:self.doc_or_operation(12))
        self.button_save.place(relx=1.0, rely=0.5, anchor=SE)
        
        self.button_clean = Button(text='Clean', command=self.clean_fields)
        self.button_clean.place(relx=1.0, rely=0.8, anchor=SE)
        
        self.button_create_chart = Button(text='Chart', command=PJU.create_plot)
        self.button_create_chart.place(relx=1.0, rely=0.65, anchor=SE)
        
        self.button_list = [self.button_randomDNASequence,self.button_sequenceIsValid,self.button_countAcidNucl,self.button_Translate_DNAtoRNA ,self.button_Translate_toProtein,    self.button_ReverseDNAComplement, self.button_Rate_GC, self.button_monoisotopicMass, self.button_codonsFrequence,self.button_RNASplicing ,self.button_DNAAssembler, self.button_save]
        
    def doc_Functions(self):
        self.help_bool = not self.help_bool
        i=0
        if self.help_bool:
            for i in range(len(self.button_list)):
                self.button_list[i].config(bg='red', fg='white')
        else:
            for i in range(len(self.button_list)):
                self.button_list[i].config(bg='white', fg='black')
            
            self.clean_fields()
    
    def doc_or_operation(self, numFun):
        
        self.numFun = numFun #Save choosed function
        
        if self.help_bool:
            self.switchDocs()
        else:
            self.switchOperations()
       
    def switchDocs(self):
        
        doc = {
            1:lambda x:PJU.randomDNASequence.__doc__,
            2:lambda x:PJU.sequenceIsValid.__doc__,
            3:lambda x:PJU.countAcidNucl.__doc__,
            4:lambda x:PJU.Translate_DNAtoRNA.__doc__,
            5:lambda x:PJU.Translate_toProtein.__doc__,
            6:lambda x:PJU.ReverseDNAComplement.__doc__,
            7:lambda x:PJU.Rate_GC.__doc__,
            8:lambda x:PJU.codonsFrequence.__doc__,
            9:lambda x:PJU.monoisotopicMass.__doc__,
            10:lambda x:PJU.RNASplicing.__doc__,
            11:lambda x:PJU.DNAAssembler.__doc__,
            12:lambda x:self.saveResult.__doc__,
            }[numFun](self.numFun)
        self.viewResult(doc)
        
    def switchOperations(self):
        """Execute the function based on numFun.
        
        This function is created to distinguich minor operation from 
        major ones. It will verify the value of numFun so if it is equal
        to 1 """
        if self.numFun == 1:
            self.createRandomSeq()
            
        elif self.numFun == 12:
            self.saveResult()
        else:
            self.open_input_window(), 
            
    def open_input_window(self):
        self.fenetre_input = Toplevel(self.master)
        self.fenetre_input.grab_set() 
        self.fenetre_input.resizable(0, 0) #To prevent resizing
        self.fenetre_input.geometry('500x140')
        self.fenetre_input.title('Input a sequence - Traitement {}'.format(self.numFun))
        
        bb = Inputwindow(self.fenetre_input, self.numFun)
        
        result = bb.show()
        self.fenetre_input.grab_release()
        self.viewResult(result)
        
    def clean_fields(self):
        self.entryText.set('')
        self.result_text.config(state='normal')
        self.result_text.delete('1.0', END)
        self.result_text.config(state='disable')
        
    def createRandomSeq(self):
        """Verify if a number was entered to create a random sequence which length equal to number.
        
        verify if Entry named 'entry_randomDNASeq' is not empty than
        call the function named 'randomDNASequence' by passing the
        value of the 'entry_randomDNASeq' as paramater to the function. """
        
        length = self.entry_randomDNASeq.get()
        result = "Please enter a number in the field in the right of the button"
        if length != '':
            result = PJU.randomDNASequence(length)
        
        self.viewResult(result)
    
    def saveResult(self):
        """Save text in result field by clicking on button labelled 
        'Save in file'.
        
        save in a file with current wd with a message in the end of the
        result. This will not erase result content. 
        Passing by a while loop to distinguich current files in saving 
        folder to not erase existant and a get a sequenced saved results
        files."""
        
        if self.result_text.get('1.0', END).strip('\n') == '':
            return
        
        if 'SUCCESSFULY SAVED' in self.result_text.get('1.0', END).strip('\n'):
            self.viewResult('\n\nALREADY SAVED', add=True)
            return
        
        try:
            while True:
                open('saving_result_{}'.format(self.file_index), 'r')
                self.file_index += 1
        except:
            filename = 'saving_result_{}.txt'.format(self.file_index)
            fi = open(filename, 'w')
            fi.write(self.result_text.get('1.0', END))
            self.viewResult('\n\nSUCCESSFULY SAVED in {}.'.format(filename), add=True)

    def viewResult(self, result, add=False):
        """View the result of an operation on the sequences by clicking on button named 'Launch' in the input window.
        
        Change the content of the widget Text named 'result_text' by
        deleting the content and replacing it with the value of 
        parameter 'text'.
        See button_Launch in file GUI_InputWindow.py."""
        self.result_text.config(state='normal')
        if not add:
            self.result_text.delete('1.0', END)
        self.result_text.insert(END, str(result))
        self.result_text.config(state='disable')
        
def ui_main_window():
    """Launch main window.
    
    """
    #main window intialize
    fenetre = Tk()
    fenetre.resizable(0, 0) #To prevent resizing
    fenetre.geometry('800x600')
    fenetre.title('Bioinformatics Application')
    interface = Interface(fenetre)
    try:
        interface.mainloop()
    except TclError:
        #To prevent error if the user close the table without clicking on button_exit
        interface.destroy()

class Inputwindow(tk.Toplevel):
    """Input sequences window.
    
    Here you will enter the file path or the DNA/RNA/Protein sequence 
    in an inputField -Entry-, you will add option thanks to 
    Radiobuttons and Checkboxes."""

    def __init__(self, master, numFun):
        self.top = self
        self.master = master
        self.frame_1 = Frame(self.master)
        self.frame_1.pack(side=BOTTOM, fill=X, pady=5)
        self.frame_2 = Frame(self.master)
        self.frame_2.pack(side=TOP, fill=BOTH, pady=5)
        self.frame_3 = Frame(self.frame_2)
        self.frame_3.pack(side=RIGHT, fill=BOTH)
        
        self.label_input = Label(self.frame_2, text='Fichier/Sequence : ')
        self.label_input.pack(side=LEFT)
        
        self.check_introns_var = IntVar()
        self.check_introns = Checkbutton(self.frame_2, text='Introns', variable=self.check_introns_var, command=self.introns_input)
        
        self.entry_sequence = Text(self.frame_2, width=40)
        self.entry_sequence.pack(side=LEFT)
        
        self.check_fromFile_var = IntVar()
        self.check_fromFile = Checkbutton(self.frame_3, text="is file ", variable=self.check_fromFile_var)
        self.check_fromFile.pack(side=TOP, pady=8)
        
        self.check_isProtein_var = IntVar()
        self.check_isProtein = Checkbutton(self.frame_3, text="is protein ", variable=self.check_isProtein_var)
        
        self.frame_4 = Frame(self.frame_3)

        self.var_choice_DNA_RNA = StringVar()
        self.choice_dna = Radiobutton(self.frame_4, text='DNA', variable=self.var_choice_DNA_RNA, value='dna')
        self.choice_rna = Radiobutton(self.frame_4, text='RNA', variable=self.var_choice_DNA_RNA, value='rna')
        self.choice_dna.pack(side=LEFT)
        self.choice_rna.pack(side=RIGHT)
        
        self.button_Launch = Button(self.frame_1, text='Launch', command=lambda:self.execute_inputs(numFun))
        self.button_Launch.pack(side=RIGHT)
        
        self.button_EXAMPLE = Button(self.frame_1, text='use a ready EXAMPLE', command=lambda:self.charge_example(numFun))
        self.button_EXAMPLE.pack(side=RIGHT, padx=100)
        
        self.button_discard = Button(self.frame_1, text='discard', command=self.master.destroy)
        self.button_discard.pack(side=LEFT)
        
#        self.label_information = Label(self.frame_2,text="Separate sequences with '>'")
#        self.label_information.place(relx=0.5, rely=0.0, anchor=CENTER)
        
        self.launch_result=''
        
        self.packSpecialWidgets(numFun)
    
    def packSpecialWidgets(self, numFun):
        
        if numFun == 10:
            self.check_introns.place(relx=0.0, rely=0.1, anchor='nw')
            self.introns_variable = ''
            self.sequence_variable = ''
            return
        if numFun in (4, 6, 7):
            return
        
        self.frame_4.pack(side=BOTTOM, fill=X, pady=5)
        
        if numFun in (2,9):
            self.check_isProtein.pack(side=TOP)       
        
       
    def introns_input(self):
        if self.check_introns_var.get() == 1:
            self.label_input.config(text='introns only :')
            self.sequence_variable = self.entry_sequence.get("1.0", "end-1c")
            self.entry_sequence.delete('1.0', END)
            self.entry_sequence.insert(END, self.introns_variable)
            self.button_EXAMPLE.config(state=DISABLED)
        else:
            self.label_input.config(text='Fichier/Sequence : ')
            self.introns_variable = self.entry_sequence.get("1.0", "end-1c")
            self.entry_sequence.delete('1.0', END)
            self.entry_sequence.insert(END, self.sequence_variable)
            self.button_EXAMPLE.config(state=NORMAL)
            
    def show(self):
        """will return the result to the main window.
        
        Main window will wait until the input window is closed than and
        only than return the result stocked in the input window variable
        'launch_result'.
        This function is only called by 'open_input_window' function in
        class 'Interface'."""
        
        self.master.wait_window()
        return self.launch_result
    
    def execute_inputs(self, numFun):
        """Execute an operation and save the result in 'self.launch_result' variable.
        
        Based on parameter 'numFun', the function will know which util 
        call in 'ProjectUtils' module which is imported as 'PJU'.
        Before begin sequence treatement, must verify if text_field is
        not empty and the variable of Radiobutton_DNA_RNA is not empty
        (choosed one), after this make a switch statement based on 
        the use of a dict and returning the result in the variable 
        'result'. Saving the result in the variable is the next and 
        last instruction.
        This function in called by clicking on button named 
        'button_Launch' labelled by 'Launch'."""
        
        seq_or_file = self.entry_sequence.get("1.0", "end-1c")
        if len(seq_or_file) == 0:#Easy way to verify if Text widget is empty
            #Show Error pop-up
            print("Show Error pop-up")
            return
        if self.check_isProtein_var.get() != 1:
            if self.var_choice_DNA_RNA.get() == '' and self.check_isProtein_var.get() == 0:
                #Show Error pop-up
                print("Show Error pop-up rna")
                return
        
        #ELSE EVerything is alright
        fromFile = True if self.check_fromFile_var.get() == 1 else False
        isProtein = True if self.check_isProtein_var.get() == 1 else False
        isDNA = True if self.var_choice_DNA_RNA.get() in 'dna' else False
        
        #This is a replacement for a switch case
        result ={
            2:lambda x:PJU.sequenceIsValid(seq_or_file, fromFile=fromFile, isProtein=isProtein, isDNA=isDNA),
            3:lambda x:self.counting_acid_nucleics(seq_or_file, isDNA, fromFile),
            4:lambda x:PJU.Translate_DNAtoRNA(seq_or_file, fromFile=fromFile),
            5:lambda x:PJU.Translate_toProtein(seq_or_file, isDNA=isDNA, fromFile=fromFile),
            6:lambda x:PJU.ReverseDNAComplement(seq_or_file, fromFile=fromFile),
            7:lambda x:PJU.Rate_GC(seq_or_file, fromFile=fromFile),
            8:lambda x:PJU.codonsFrequence(seq_or_file, isDNA=isDNA, fromFile=fromFile),
            9:lambda x:PJU.monoisotopicMass(seq_or_file, isDNA=isDNA, fromFile=fromFile, isProtein=isProtein),
            10:lambda x:PJU.RNASplicing(seq_or_file, self.introns_variable, fromFile=fromFile),
            11:lambda x:PJU.DNAAssembler(seq_or_file, fromFile=fromFile)
        }[numFun](numFun)
        self.launch_result = result
        self.master.destroy()
        
    def counting_acid_nucleics(self, seq, isDNA, fromFile):
        """Manage inputs to add count a specific acid Nucleic.
        
        This function is called by function named 'execute_inputs' by
        choosing numFun equal to 3 after clicking on button in the 
        Main window called 'button_countAcidNucl'.
        """
        return PJU.countAcidNucl(seq, acid=None, isDNA=isDNA, fromFile=fromFile)
   
    def charge_example(self, numFun):
        """Charge prepared example.
        
        If too lazy to enter sequences which is the case with all good 
        programmers, clic here.
        This will generate ready DNA sequences to use. Make your self
        confortable."""
        
        if numFun == 10:
            self.introns_variable = '>Intron1\nATCGGTCGAA\n>Intron2\nATCGGTCGAGCGTGT'
        
        if numFun == 11:
            exp = '>Sequence1\nATTAGACCTG\n>Sequence2\nCCTGCCGGAA\n>Sequence3\nAGACCTGCCG\n>Sequence4\nGCCGGAATAC'
        else:
            exp =[ '>Rosalind_56\nATTAGACCTG\n>Rosalind_57\nCCTGCCGGAA\n>Rosalind_58\nAGACCTGCCG\n>Rosalind_59\nGCCGGAATAC', 'ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG', PJU.randomDNASequence(random.randint(30, 300))]
            
            c = random.randint(0,len(exp)-1)
            exp = exp[c]
        
        self.entry_sequence.delete('1.0', END)
        self.entry_sequence.insert(END, exp)
        self.check_fromFile_var.set(0)
        self.check_isProtein_var.set(0)
        self.var_choice_DNA_RNA.set('dna')
                
        return True

if __name__ == "__main__":
    
    ui_main_window()
