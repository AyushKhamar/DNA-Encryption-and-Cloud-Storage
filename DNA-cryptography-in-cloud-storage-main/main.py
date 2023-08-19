from tkinter import *
from tkinter import filedialog
import os
import requests
from dna_get import main as encrypt
from dna_gdt import main as decrypt
import requests
import json
from tkinter import messagebox
from drive import upload_file, download_file, search_file
root = Tk()
root.geometry("500x500")
root.title(
    "DNA Based Encryption and decryption of confidential data in cloud storage ! ")
temp_path = r'C:\Users\DELL\OneDrive\Desktop\DNA-cryptography-in-cloud-storage-main (2)\DNA-cryptography-in-cloud-storage-main\DNA-cryptography-in-cloud-storage-main\encrypted.txt'


class encryptor:
    def __init__(self):
        self.path = None
        self.folder = None

    def browsefile(self):
        filename = filedialog.askopenfilename()
        return filename

    def getfilepath(self):
        self.path = self.browsefile()
        return self.path


def send_data_to_cloud(text):
    x = requests.get('http://localhost/phpfiles/store.php?text='+text)
    if x.text == "Success":
        messagebox.showinfo(
            "Success", "Encryption Successfull!! Please enter Google Account for google drive upload")

    else:
        messagebox.showerror("Error", "Pls check weather the host is up !")


def get_data_from_cloud():
    y = requests.get('http://localhost/phpfiles/retrieve.php')
    x = json.loads(y.text)
    return x['data']


def Encryption_upload():
    q = encryptor()
    textdata = q.browsefile()
    encrypted_string = encrypt(textdata)
    send_data_to_cloud(encrypted_string)
    messagebox.showinfo("Encrypted text", encrypted_string)
    messagebox.showwarning(
        "KEY", "A key.txt will be generated, please keep it safe because you will need it for decryption ")
    upload_file('encrypted.txt', temp_path, 'text/plain')


def download_decryption():

    encrypted_text = get_data_from_cloud()
    decrypted_string = decrypt(encrypted_text)
    messagebox.showinfo("Decrypted data", decrypted_string)
    messagebox.showinfo(
        "Decrypted data", "Decryption Successfull!! , please check decrypted.txt to get the decrypted information")


l1 = Label(
    root, text="Pls select any file to encrypt and upload to the cloud ! ", fg="green")

upload_button = Button(text="Encrypt and Upload a File", height=3,
                   width=50, command=Encryption_upload, fg="white", bg="green")

download_button = Button(text="Download and decrypt the file", height=3,
                     width=50, command=download_decryption, fg="white", bg="black")
l1.pack()
upload_button.pack()
download_button.pack()
l1.place(relx=0.5, rely=0.3, anchor=CENTER)
upload_button.place(relx=0.5, rely=0.4, anchor=CENTER)
download_button.place(relx=0.5, rely=0.5, anchor=CENTER)


root.mainloop()
