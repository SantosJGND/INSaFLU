# -*- coding: utf-8 -*-
# Generated by Django 1.11.29 on 2022-12-09 10:56
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('pathogen_identification', '0008_auto_20221203_1802'),
        ('settings', '0013_parameter_dataset'),
    ]

    operations = [
        migrations.AddField(
            model_name='parameter',
            name='televir_project',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.PROTECT, related_name='parameter', to='pathogen_identification.Projects'),
        ),
        migrations.AlterField(
            model_name='parameter',
            name='name',
            field=models.CharField(blank=True, db_index=True, max_length=100, null=True),
        ),
        migrations.AlterField(
            model_name='parameter',
            name='parameter',
            field=models.CharField(blank=True, db_index=True, max_length=150, null=True),
        ),
    ]