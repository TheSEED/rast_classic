#
# Copyright (c) 2003-2006 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
# 
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License. 
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#

package DBrtns;

    # Inherit the DBKernel methods. We must do this BEFORE the "use strict".
    use DBKernel;
    @ISA = qw(DBKernel);

use strict;
use POSIX;
use DBI;
use FIG_Config;

use Data::Dumper;
use Carp;

sub new {
    my($class,$dbms,$dbname,$dbuser,$dbpass,$dbport, $dbhost, $dbsock) = @_;

    $dbms   = defined($dbms)   ? $dbms   : $FIG_Config::dbms;
    $dbname = defined($dbname) ? $dbname : $FIG_Config::db;
    $dbuser = defined($dbuser) ? $dbuser : $FIG_Config::dbuser;
    $dbpass = defined($dbpass) ? $dbpass : $FIG_Config::dbpass;
    $dbport = defined($dbport) ? $dbport : $FIG_Config::dbport;
    $dbhost = defined($dbhost) ? $dbhost : $FIG_Config::dbhost;
    $dbsock = defined($dbsock) ? $dbsock : $FIG_Config::dbsock;

    return DBKernel::new($class, $dbms, $dbname, $dbuser, $dbpass, $dbport, $dbhost, $dbsock);
}

=head1 get_inserted_id

Return the last ID of a row inserted into an autonumber/serial-containing table.

=cut

sub get_inserted_id {
    my($self, $table, $sth, $id_column) = @_;

    $id_column = 'id' unless defined($id_column);
    if ($self->{_dbms} eq "Pg") {
        my $oid = $sth->{pg_oid_status};
        my $ret = $self->SQL("select $id_column from $table where oid = ?", undef, $oid);
        return $ret->[0]->[0];
    } elsif ($self->{_dbms} eq "mysql") {
        my $id = $self->{_dbh}->{mysql_insertid};
        # print "mysql got $id\n";
        return $id;
    }
    else
    {
	confess "Attempting get_inserted_id on unsupported database $self->{_dbms}\n";
    }
}

#
# Following are database administration routines. They create an instance of a ServerAdmin class
# for the appropriate server type (in order to eliminate the if mysql / if pg / etc stuff).
#

sub get_server_admin
{
    if ($FIG_Config::dbms eq "mysql")
    {
	return MysqlAdmin->new();
    }
    elsif ($FIG_Config::dbms eq "Pg")
    {
	return new PostgresAdmin();
    }
    else
    {
	warn "Unknown server type $FIG_Config::dbms\n";
	return undef;
    }
}
package MysqlAdmin;

use POSIX;
use DBI;

sub new
{
    my($class) = @_;

    my $self = {};

    return bless($self, $class);
}

sub init_db
{
    my($self, $db_dir) = @_;

    if (!$db_dir)
    {
	warn "init_db failed: db_dir must be provided\n";
	return;
    }

    if (-d "$db_dir/mysql")
    {
	warn "init_db: mysql data directory already exists\n";
	return;
    }

    my $exe = "$FIG_Config::ext_bin/mysql_install_db";
    if (! -x $exe)
    {
	$exe = "mysql_install_db";
    }


    my @opts;

    push(@opts, "--datadir=$db_dir");
    push(@opts, "--user=$FIG_Config::dbuser");

    if (not $FIG_Config::use_system_mysql)
    {
	push(@opts, "--basedir=$FIG_Config::common_runtime")
    }


    my $rc = system($exe, @opts);
    if ($rc != 0)
    {
	my $err = $?;
	if (WIFEXITED($err))
	{
	    my $exitstat = WEXITSTATUS($err);
	    warn "init_db failed: $exe returned result code $exitstat\n";
	}
	else
	{
	    warn "init_db failed: $exe died with signal ", WTERMSIG($err), "\n";
	}
	return;
    }

    return 1;
}

sub create_database
{
    my($self, $db_name) = @_;

    my $drh = DBI->install_driver("mysql");

    my @dbs = DBI->data_sources("mysql", { host => $FIG_Config::dbhost,
					       user => $FIG_Config::dbuser,
					       password => $FIG_Config::dbpass });
    if (grep { $_ eq $db_name } @dbs)
    {
	warn "Database $db_name already exists\n";
	return;
    }

    my $rc = $drh->func('createdb', $db_name, $FIG_Config::dbhost,
			$FIG_Config::dbuser, $FIG_Config::dbpass, 'admin');


    if (!$rc)
    {
	warn "create_database: createdb call failed: $DBI::errstr\n";
	return;
    }

    return 1;
}

sub start_server
{
    my($self, $dont_fork) = @_;

    print "Starting mysql server\n";

    my(@opts);

    my $cnf = "$FIG_Config::fig_disk/config/my.cnf";

    if ($FIG_Config::use_system_mysql)
    {
	#
	# This has to be first in the argument list.
	#
	push(@opts, "--defaults-extra-file=$cnf");
    }

    #
    # Put this first, so  config can put --defaults-extra-file here 
    # and have it show up first.
    #
    if (@FIG_Config::db_server_startup_options)
    {
	push(@opts, @FIG_Config::db_server_startup_options)
    }

    push(@opts, "--port=$FIG_Config::dbport");
    #
    # Don't do this; dbuser isn't the unix uid that we are using.
    #
    #push(@opts, "--user=$FIG_Config::dbuser");

    push(@opts, "--datadir=$FIG_Config::db_datadir");

    if ($FIG_Config::use_system_mysql)
    {
	push(@opts, "--err-log=$FIG_Config::temp/mysql.log");
	push(@opts, "--socket=$FIG_Config::dbsock");

	#
	# Feh. You can't actually override the socket that /etc/my.cnf 
	# sets up, so we need to set up a config/my.cnf with the socket in it.
	#

	if (! -f $cnf) 
	{
	    if (open(F, ">$cnf"))
	    {
		print F <<END;
[mysqld]
socket=$FIG_Config::dbsock
END
		close(F);
	    }
	}
    }
    else
    {
	push(@opts, "--basedir=$FIG_Config::common_runtime");
	push(@opts, "--ledir=$FIG_Config::common_runtime/libexec");
    }

    if (not $FIG_Config::mysql_v3)
    {
	push(@opts, "--old-password");
	push(@opts, "--max-allowed-packet=128M");
    }

    #
    # Use InnoDB for large-table support and allegedly better performance.
    #
    
    #push(@opts, "--default-table-type=innodb");
    
    #
    # Oddly, this doesn't seem to work. need to set the environment variable.
    #
    #push(@opts, "--port=$FIG_Config::dbport");

    #
    # We are going to assume that if mysql has shipped with this release, we'll use it. Otherwise
    # try to use a system one.
    #

    my $exe;
    if ($FIG_Config::mysql_v3)
    {
	$exe = "safe_mysqld";
    }
    else
    {
	$exe = "mysqld_safe";
    }

    if (-x "$FIG_Config::ext_bin/$exe")
    {
	$exe = "$FIG_Config::ext_bin/$exe";
    }

    print "Start $exe @opts\n";

    if ($dont_fork)
    {
	$ENV{MYSQL_TCP_PORT} = $FIG_Config::dbport;
	exec $exe, @opts;
    }
    else
    {
	my $pid = fork;

	if ($pid == 0)
	{
	    POSIX::setsid();

	    $ENV{MYSQL_TCP_PORT} = $FIG_Config::dbport;
	    exec $exe, @opts;
	}
	print "Forked db server $pid\n";
    }

}

1;
