Deployment Script
=================

*INSTALLER* : The machine from where the deployment script will be run

*SERVER*    : The main server where we want to deploy the project

This setup requires `ansible` to be installed in the INSTALLER.
If the machine has OSX as the operating system and `homebrew` is installed,
please use this command:

```
$ brew install ansible
```

Add an entry in the `/etc/ansible/hosts` file (create it if not preasent)
to give detail of the deployment server:

```
vdc4ml ansible_host=100.100.100.100 ansible_connection=ssh ansible_user=root
# here 'vdc4ml' is the given name of the machine
# IP address 100.100.100.100
# Connect using ssh
# Will be connecting as root
```

The RSA key-pair need to be generated on *INSTALLER* and the public key needs
to be copied to the `~/.ssh/authorized_keys` file on the *SERVER*.
