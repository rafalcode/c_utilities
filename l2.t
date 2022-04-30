Goal
By the end of this tutorial, you should have a JupyterHub with some admin users and a user environment with packages you want installed running on a server you have access to.

Pre-requisites
Some familiarity with the command line.

A server running Ubuntu 18.04 where you have root access.

At least 1GB of RAM on your server.

Ability to ssh into the server & run commands from the prompt.

An IP address where the server can be reached from the browsers of your target audience.

If you run into issues, look at the specific troubleshooting guide for custom server installations.

Step 1: Installing The Littlest JupyterHub
Using a terminal program, SSH into your server. This should give you a prompt where you can type commands.

Make sure you have python3, python3-dev, curl and git installed.

sudo apt install python3 python3-dev git curl
￼
Copy the text below, and paste it into the terminal. Replace <admin-user-name> with the name of the first admin user for this JupyterHub. Choose any name you like (don’t forget to remove the brackets!). This admin user can log in after the JupyterHub is set up, and can configure it to their needs. Remember to add your username!

curl -L https://tljh.jupyter.org/bootstrap.py | sudo -E python3 - --admin <admin-user-name>
